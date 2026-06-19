//! Evidence clusterer — port of the Scala `SvEvidenceClusterer`. Clusters discordant
//! pairs + split reads into SV calls: inter-chromosomal -> BND, intra-chromosomal
//! positional clustering with orientation/insert-based type inference, then integration
//! with depth-based CNV segments.

use std::collections::{BTreeMap, BTreeSet};

use super::evidence::{
    BreakpointCluster, DepthSegment, DiscordantPair, DiscordantReason, SplitRead, SvEvidenceCollection,
};
use super::segmenter;
use super::types::{SvCall, SvCallerConfig, SvType};

enum EvidencePoint {
    Pair(DiscordantPair),
    Split(SplitRead),
}

/// Cluster all SV evidence into calls, integrating depth segments.
pub fn cluster(
    evidence: &SvEvidenceCollection,
    depth_segments: &[DepthSegment],
    config: &SvCallerConfig,
) -> Vec<SvCall> {
    let mut calls: Vec<SvCall> = Vec::new();
    let mut call_index = 0;

    // 1. Translocations (inter-chromosomal).
    for cluster in cluster_translocations(&evidence.inter_chromosomal_pairs(), config) {
        call_index += 1;
        calls.push(breakpoint_cluster_to_call(&cluster, SvType::Bnd, call_index, config));
    }

    // 2. Intra-chromosomal events.
    let intra: Vec<&DiscordantPair> = evidence
        .discordant_pairs
        .iter()
        .filter(|p| p.reason != DiscordantReason::InterChromosomal)
        .collect();

    let mut pairs_by_chrom: BTreeMap<&str, Vec<DiscordantPair>> = BTreeMap::new();
    for p in &intra {
        pairs_by_chrom.entry(p.chrom1.as_str()).or_default().push((*p).clone());
    }
    let mut splits_by_chrom: BTreeMap<&str, Vec<SplitRead>> = BTreeMap::new();
    for s in &evidence.split_reads {
        splits_by_chrom
            .entry(s.primary_chrom.as_str())
            .or_default()
            .push(s.clone());
    }
    let chroms: BTreeSet<&str> = pairs_by_chrom.keys().chain(splits_by_chrom.keys()).copied().collect();

    for chrom in chroms {
        let pairs = pairs_by_chrom.get(chrom).cloned().unwrap_or_default();
        let splits = splits_by_chrom.get(chrom).cloned().unwrap_or_default();
        for cluster in cluster_intra(chrom, pairs, splits, config) {
            if cluster.total_support() >= config.min_total_support {
                let sv_type = infer_sv_type(&cluster);
                call_index += 1;
                calls.push(breakpoint_cluster_to_call(&cluster, sv_type, call_index, config));
            }
        }
    }

    // 3. Integrate with depth segments, then sort.
    let mut integrated = integrate_pe_sr_with_depth(calls, depth_segments, config);
    integrated.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
    integrated
}

/// Cluster inter-chromosomal pairs by chromosome-pair, then by position on chrom1.
fn cluster_translocations(pairs: &[DiscordantPair], config: &SvCallerConfig) -> Vec<BreakpointCluster> {
    let mut by_chrom_pair: BTreeMap<(&str, &str), Vec<DiscordantPair>> = BTreeMap::new();
    for p in pairs {
        by_chrom_pair
            .entry((p.chrom1.as_str(), p.chrom2.as_str()))
            .or_default()
            .push(p.clone());
    }

    let mut out = Vec::new();
    for ((chrom1, chrom2), group) in by_chrom_pair {
        let positioned: Vec<(i64, DiscordantPair)> = group.iter().map(|p| (p.pos1, p.clone())).collect();
        for (pos, clustered) in cluster_by_position(positioned, config) {
            let mate_positions: Vec<i64> = clustered.iter().map(|(_, p)| p.pos2).collect();
            let mate_pos = mate_positions.iter().sum::<i64>() / mate_positions.len() as i64;
            out.push(BreakpointCluster {
                chrom: chrom1.to_string(),
                position: pos,
                ci_low: -(config.max_cluster_distance as i32),
                ci_high: config.max_cluster_distance as i32,
                discordant_pairs: clustered.into_iter().map(|(_, p)| p).collect(),
                split_reads: Vec::new(),
                mate_chrom: Some(chrom2.to_string()),
                mate_position: Some(mate_pos),
            });
        }
    }
    out
}

/// Cluster (position, item) pairs greedily within `max_cluster_distance` of the cluster start.
fn cluster_by_position(
    items: Vec<(i64, DiscordantPair)>,
    config: &SvCallerConfig,
) -> Vec<(i64, Vec<(i64, DiscordantPair)>)> {
    if items.is_empty() {
        return Vec::new();
    }
    let mut sorted = items;
    sorted.sort_by_key(|(p, _)| *p);

    let mut clusters = Vec::new();
    let mut current: Vec<(i64, DiscordantPair)> = Vec::new();
    let mut cluster_start = sorted[0].0;
    for item in sorted {
        if item.0 - cluster_start <= config.max_cluster_distance || current.is_empty() {
            current.push(item);
        } else {
            let mean = current.iter().map(|(p, _)| *p).sum::<i64>() / current.len() as i64;
            clusters.push((mean, std::mem::take(&mut current)));
            cluster_start = item.0;
            current.push(item);
        }
    }
    if !current.is_empty() {
        let mean = current.iter().map(|(p, _)| *p).sum::<i64>() / current.len() as i64;
        clusters.push((mean, current));
    }
    clusters
}

/// Cluster intra-chromosomal PE+SR evidence by position.
fn cluster_intra(
    chrom: &str,
    pairs: Vec<DiscordantPair>,
    splits: Vec<SplitRead>,
    config: &SvCallerConfig,
) -> Vec<BreakpointCluster> {
    let mut points: Vec<(i64, EvidencePoint)> = Vec::new();
    for p in pairs {
        points.push((p.pos1, EvidencePoint::Pair(p)));
    }
    for s in splits {
        points.push((s.primary_pos, EvidencePoint::Split(s)));
    }
    if points.is_empty() {
        return Vec::new();
    }
    points.sort_by_key(|(p, _)| *p);

    let mut clusters = Vec::new();
    let mut current: Vec<(i64, EvidencePoint)> = Vec::new();
    let mut cluster_start = points[0].0;
    for point in points {
        if point.0 - cluster_start <= config.max_cluster_distance || current.is_empty() {
            current.push(point);
        } else {
            clusters.push(make_cluster(chrom, std::mem::take(&mut current)));
            cluster_start = point.0;
            current.push(point);
        }
    }
    if !current.is_empty() {
        clusters.push(make_cluster(chrom, current));
    }
    clusters
}

fn make_cluster(chrom: &str, points: Vec<(i64, EvidencePoint)>) -> BreakpointCluster {
    let positions: Vec<i64> = points.iter().map(|(p, _)| *p).collect();
    let mean_pos = positions.iter().sum::<i64>() / positions.len() as i64;
    let min_pos = *positions.iter().min().unwrap();
    let max_pos = *positions.iter().max().unwrap();

    let mut pairs = Vec::new();
    let mut splits = Vec::new();
    for (_, ev) in points {
        match ev {
            EvidencePoint::Pair(p) => pairs.push(p),
            EvidencePoint::Split(s) => splits.push(s),
        }
    }
    BreakpointCluster {
        chrom: chrom.to_string(),
        position: mean_pos,
        ci_low: (min_pos - mean_pos) as i32,
        ci_high: (max_pos - mean_pos) as i32,
        discordant_pairs: pairs,
        split_reads: splits,
        mate_chrom: None,
        mate_position: None,
    }
}

/// Infer SV type from a cluster's evidence pattern.
fn infer_sv_type(cluster: &BreakpointCluster) -> SvType {
    if cluster.mate_chrom.is_some() {
        return SvType::Bnd;
    }
    let n = cluster.discordant_pairs.len();
    let same_strand = cluster
        .discordant_pairs
        .iter()
        .filter(|p| p.strand1 == p.strand2)
        .count();
    if same_strand > n / 2 {
        return SvType::Inv;
    }
    let outliers: Vec<&DiscordantPair> = cluster
        .discordant_pairs
        .iter()
        .filter(|p| p.reason == DiscordantReason::InsertSizeOutlier)
        .collect();
    if !outliers.is_empty() {
        let avg_insert = outliers.iter().map(|p| p.insert_size as f64).sum::<f64>() / outliers.len() as f64;
        let expected_insert = cluster
            .discordant_pairs
            .first()
            .map_or(400.0, |p| p.insert_size as f64 * 0.5);
        if avg_insert > expected_insert * 2.0 {
            return SvType::Del;
        } else if avg_insert < expected_insert * 0.5 {
            return SvType::Dup;
        }
    }
    SvType::Del
}

fn breakpoint_cluster_to_call(
    cluster: &BreakpointCluster,
    sv_type: SvType,
    index: i32,
    config: &SvCallerConfig,
) -> SvCall {
    let quality = (cluster.total_support() as f64 * 5.0 + cluster.mean_mapq() * 0.5).min(99.0);

    let (sv_len, end, mate_chrom, mate_pos) = if sv_type == SvType::Bnd {
        (
            0i64,
            cluster.position,
            cluster.mate_chrom.clone(),
            cluster.mate_position,
        )
    } else {
        let mate_positions: Vec<i64> = cluster
            .discordant_pairs
            .iter()
            .map(|p| p.pos2)
            .filter(|&p| p != cluster.position)
            .collect();
        if !mate_positions.is_empty() {
            let avg = mate_positions.iter().sum::<i64>() / mate_positions.len() as i64;
            let length = (avg - cluster.position).abs();
            let signed = if sv_type == SvType::Del { -length } else { length };
            (signed, cluster.position + length, None, None)
        } else {
            (
                config.max_cluster_distance,
                cluster.position + config.max_cluster_distance,
                None,
                None,
            )
        }
    };

    let genotype = if cluster.total_support() >= 10 { "1/1" } else { "0/1" };
    let filter = if cluster.pe_support() >= config.min_paired_end_support
        || cluster.sr_support() >= config.min_split_read_support
    {
        "PASS"
    } else {
        "LowSupport"
    };

    SvCall {
        id: format!("{}_{}_{}_{}", sv_type.as_str(), cluster.chrom, cluster.position, index),
        chrom: cluster.chrom.clone(),
        start: cluster.position,
        end,
        sv_type,
        sv_len,
        ci_pos: (cluster.ci_low, cluster.ci_high),
        ci_end: (cluster.ci_low, cluster.ci_high),
        quality,
        paired_end_support: cluster.pe_support(),
        split_read_support: cluster.sr_support(),
        relative_depth: None,
        mate_chrom,
        mate_pos,
        filter: filter.into(),
        genotype: genotype.into(),
    }
}

/// Attach depth evidence to overlapping PE/SR calls; append depth-only calls.
fn integrate_pe_sr_with_depth(
    pe_sr_calls: Vec<SvCall>,
    depth_segments: &[DepthSegment],
    config: &SvCallerConfig,
) -> Vec<SvCall> {
    if depth_segments.is_empty() {
        return pe_sr_calls;
    }
    let depth_calls = segmenter::to_sv_calls(depth_segments, config);
    let mut used = vec![false; depth_calls.len()];

    let mut enhanced: Vec<SvCall> = pe_sr_calls
        .into_iter()
        .map(|mut call| {
            if let Some(idx) = depth_calls.iter().enumerate().position(|(idx, dc)| {
                !used[idx]
                    && dc.chrom == call.chrom
                    && dc.sv_type == call.sv_type
                    && overlaps(call.start, call.end, dc.start, dc.end)
            }) {
                used[idx] = true;
                call.relative_depth = depth_calls[idx].relative_depth;
            }
            call
        })
        .collect();

    for (idx, dc) in depth_calls.into_iter().enumerate() {
        if !used[idx] {
            enhanced.push(dc);
        }
    }
    enhanced
}

fn overlaps(start1: i64, end1: i64, start2: i64, end2: i64) -> bool {
    start1 <= end2 && start2 <= end1
}
