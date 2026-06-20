//! Y-STR autoclustering with SNP-branch propagation (FTDNA project-import follow-on).
//!
//! Given a project's members — each with a Y-STR haplotype and, for the SNP-placed ones, a branch
//! label — group them into clusters by Y-STR genetic distance, then **propagate** each cluster's
//! SNP branch onto its STR-only members as a *suggested* placement. Effectively: "this STR-only
//! haplotype most likely sits on SNP branch X."
//!
//! Pure marker math, no IO. Genetic distance is order-independent for multi-copy markers and
//! normalized per 100 comparable markers so mixed panel sizes (Y-12 … Y-700) compare fairly.

use std::collections::HashMap;

use du_domain::ids::SampleGuid;
use serde::{Deserialize, Serialize};

use crate::strprofile::StrMarker;

/// One project member fed to the clusterer.
#[derive(Debug, Clone)]
pub struct ClusterMember {
    pub guid: SampleGuid,
    /// Display label (kit / name).
    pub label: String,
    /// Confirmed SNP branch (terminal label), or `None` for an STR-only member.
    pub branch: Option<String>,
    pub markers: Vec<StrMarker>,
}

/// Tuning for [`cluster_ystr`].
#[derive(Debug, Clone)]
pub struct ClusterOpts {
    /// Minimum shared markers for two members to be comparable.
    pub min_markers: i64,
    /// Max normalized genetic distance (mutations per 100 markers) to link two members.
    pub link_gd_per_100: f32,
}

impl Default for ClusterOpts {
    fn default() -> Self {
        // ~8 mutations / 100 markers groups close kin/sub-branches without chaining the whole clade.
        Self {
            min_markers: 25,
            link_gd_per_100: 8.0,
        }
    }
}

/// A suggested branch for an STR-only member, from its nearest SNP-placed cluster-mate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BranchSuggestion {
    pub branch: String,
    /// Genetic distance (differing markers) to the nearest placed member.
    pub gd: i64,
    /// Markers compared with that nearest placed member.
    pub compared: i64,
    /// 0..1 — higher = closer match + more agreement among the cluster's placed members.
    pub confidence: f32,
}

/// A member as placed in the clustering output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClusteredMember {
    pub guid: SampleGuid,
    pub label: String,
    /// Confirmed branch (None for STR-only).
    pub branch: Option<String>,
    /// Suggested branch (only for STR-only members that landed near a placed member).
    pub suggested: Option<BranchSuggestion>,
    pub markers: usize,
}

impl ClusteredMember {
    /// The branch to display: confirmed, else suggested.
    pub fn effective_branch(&self) -> Option<&str> {
        self.branch
            .as_deref()
            .or(self.suggested.as_ref().map(|s| s.branch.as_str()))
    }
}

/// One Y-STR cluster.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct YstrCluster {
    /// Representative branch (most common confirmed branch among members), if any.
    pub branch: Option<String>,
    /// Members, placed (confirmed) first then suggested then unplaced, each sorted by label.
    pub members: Vec<ClusteredMember>,
}

impl YstrCluster {
    pub fn confirmed_count(&self) -> usize {
        self.members.iter().filter(|m| m.branch.is_some()).count()
    }
    pub fn suggested_count(&self) -> usize {
        self.members
            .iter()
            .filter(|m| m.branch.is_none() && m.suggested.is_some())
            .count()
    }
}

/// The clustering result.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct YstrClustering {
    /// Clusters, largest first.
    pub clusters: Vec<YstrCluster>,
    /// Members with too few markers to cluster (shown separately).
    pub unclustered: Vec<ClusteredMember>,
}

/// Union-find for single-linkage clustering.
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}
impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        // Path compression.
        let mut c = x;
        while self.parent[c] != r {
            let next = self.parent[c];
            self.parent[c] = r;
            c = next;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb,
            std::cmp::Ordering::Greater => self.parent[rb] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
    }
}

/// Order-independent normalized value for a (possibly multi-copy) marker (`"15-14"` ≡ `"14-15"`).
fn norm_value(v: &str) -> String {
    let mut parts: Vec<&str> = v.split('-').map(|p| p.trim()).filter(|p| !p.is_empty()).collect();
    parts.sort_unstable();
    parts.join("-")
}

/// Marker name → normalized value map (uppercased names) for fast pairwise distance.
fn marker_map(markers: &[StrMarker]) -> HashMap<String, String> {
    markers
        .iter()
        .map(|m| (m.marker.to_ascii_uppercase(), norm_value(&m.value)))
        .collect()
}

/// `(differing, compared)` over shared markers between two precomputed maps.
fn map_distance(a: &HashMap<String, String>, b: &HashMap<String, String>) -> (i64, i64) {
    let (small, large) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let (mut diff, mut comp) = (0i64, 0i64);
    for (k, va) in small {
        if let Some(vb) = large.get(k) {
            comp += 1;
            if va != vb {
                diff += 1;
            }
        }
    }
    (diff, comp)
}

/// Cluster members by Y-STR genetic distance and propagate SNP branches to STR-only members.
pub fn cluster_ystr(members: &[ClusterMember], opts: &ClusterOpts) -> YstrClustering {
    let n = members.len();
    let maps: Vec<HashMap<String, String>> = members.iter().map(|m| marker_map(&m.markers)).collect();
    // Comparable = enough markers to take part in clustering.
    let comparable: Vec<bool> = members
        .iter()
        .map(|m| (m.markers.len() as i64) >= opts.min_markers)
        .collect();

    // Single-linkage union-find over pairs within the normalized-GD threshold.
    let mut uf = UnionFind::new(n);
    for i in 0..n {
        if !comparable[i] {
            continue;
        }
        for j in (i + 1)..n {
            if !comparable[j] {
                continue;
            }
            let (diff, comp) = map_distance(&maps[i], &maps[j]);
            if comp >= opts.min_markers {
                let scaled = diff as f32 * 100.0 / comp as f32;
                if scaled <= opts.link_gd_per_100 {
                    uf.union(i, j);
                }
            }
        }
    }

    // Group comparable members by cluster root; non-comparable → unclustered.
    let mut by_root: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut unclustered = Vec::new();
    for (i, &is_comparable) in comparable.iter().enumerate() {
        if is_comparable {
            let r = uf.find(i);
            by_root.entry(r).or_default().push(i);
        } else {
            unclustered.push(clustered(members, i, None));
        }
    }

    let mut clusters: Vec<YstrCluster> = by_root
        .into_values()
        .map(|idxs| build_cluster(members, &maps, &idxs))
        .collect();
    // Largest clusters first; stable on branch label for determinism.
    clusters.sort_by(|a, b| {
        b.members
            .len()
            .cmp(&a.members.len())
            .then_with(|| a.branch.cmp(&b.branch))
    });
    unclustered.sort_by(|a, b| a.label.cmp(&b.label));

    YstrClustering { clusters, unclustered }
}

/// Assemble one cluster: propagate branches to STR-only members from their nearest placed cluster-mate.
fn build_cluster(members: &[ClusterMember], maps: &[HashMap<String, String>], idxs: &[usize]) -> YstrCluster {
    // Placed (confirmed-branch) members in this cluster.
    let placed: Vec<usize> = idxs.iter().copied().filter(|&i| members[i].branch.is_some()).collect();
    // Branch agreement: how many placed members carry each branch.
    let mut branch_counts: HashMap<&str, usize> = HashMap::new();
    for &p in &placed {
        if let Some(b) = members[p].branch.as_deref() {
            *branch_counts.entry(b).or_default() += 1;
        }
    }
    let total_placed = placed.len();

    let mut out: Vec<ClusteredMember> = idxs
        .iter()
        .map(|&i| {
            if members[i].branch.is_some() {
                return clustered(members, i, None);
            }
            // STR-only: suggest the branch of the nearest placed cluster-mate.
            let nearest = placed
                .iter()
                .filter_map(|&p| {
                    let (diff, comp) = map_distance(&maps[i], &maps[p]);
                    members[p].branch.as_deref().map(|b| (b, diff, comp))
                })
                .min_by(|a, b| a.1.cmp(&b.1).then(b.2.cmp(&a.2)));
            let suggestion = nearest.map(|(branch, gd, comp)| {
                // Confidence: closeness (GD 0 → 1.0, decaying) × branch agreement in the cluster.
                let closeness = (1.0 - gd as f32 / 12.0).clamp(0.0, 1.0);
                let agreement = if total_placed > 0 {
                    branch_counts.get(branch).copied().unwrap_or(0) as f32 / total_placed as f32
                } else {
                    0.0
                };
                BranchSuggestion {
                    branch: branch.to_string(),
                    gd,
                    compared: comp,
                    confidence: (0.5 + 0.5 * closeness) * (0.5 + 0.5 * agreement),
                }
            });
            clustered(members, i, suggestion)
        })
        .collect();

    // Confirmed first, then suggested, then unplaced; each alphabetical.
    out.sort_by(|a, b| rank(a).cmp(&rank(b)).then_with(|| a.label.cmp(&b.label)));

    // Representative branch = most common confirmed branch.
    let branch = branch_counts
        .into_iter()
        .max_by_key(|(_, c)| *c)
        .map(|(b, _)| b.to_string());

    YstrCluster { branch, members: out }
}

fn rank(m: &ClusteredMember) -> u8 {
    if m.branch.is_some() {
        0
    } else if m.suggested.is_some() {
        1
    } else {
        2
    }
}

fn clustered(members: &[ClusterMember], i: usize, suggested: Option<BranchSuggestion>) -> ClusteredMember {
    ClusteredMember {
        guid: members[i].guid,
        label: members[i].label.clone(),
        branch: members[i].branch.clone(),
        suggested,
        markers: members[i].markers.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use uuid::Uuid;

    fn m(label: &str, branch: Option<&str>, markers: &[(&str, &str)]) -> ClusterMember {
        ClusterMember {
            guid: SampleGuid(Uuid::new_v4()),
            label: label.into(),
            branch: branch.map(|s| s.to_string()),
            markers: markers
                .iter()
                .map(|(k, v)| StrMarker {
                    marker: (*k).into(),
                    value: (*v).into(),
                })
                .collect(),
        }
    }

    /// A 30-marker base haplotype; `tweak` mutates `n` markers to create genetic distance.
    fn haplo(tweak: usize) -> Vec<(&'static str, &'static str)> {
        let base = [
            ("DYS393", "13"),
            ("DYS390", "24"),
            ("DYS19", "14"),
            ("DYS391", "10"),
            ("DYS385", "11-14"),
            ("DYS426", "12"),
            ("DYS388", "12"),
            ("DYS439", "11"),
            ("DYS389I", "13"),
            ("DYS392", "13"),
            ("DYS389II", "29"),
            ("DYS458", "17"),
            ("DYS459", "9-10"),
            ("DYS455", "11"),
            ("DYS454", "11"),
            ("DYS447", "24"),
            ("DYS437", "15"),
            ("DYS448", "19"),
            ("DYS449", "29"),
            ("DYS464", "15-15-17-17"),
            ("DYS460", "11"),
            ("YGATAH4", "11"),
            ("YCAII", "19-23"),
            ("DYS456", "15"),
            ("DYS607", "15"),
            ("DYS576", "19"),
            ("DYS570", "17"),
            ("CDY", "36-38"),
            ("DYS442", "12"),
            ("DYS438", "12"),
        ];
        // Bump the value of the first `tweak` markers by appending '9' to force a difference.
        let mutated = ["99", "98", "97", "96", "95", "94", "93", "92", "91", "90", "89", "88"];
        base.iter()
            .enumerate()
            .map(|(idx, (k, v))| if idx < tweak { (*k, mutated[idx]) } else { (*k, *v) })
            .collect()
    }

    #[test]
    fn propagates_branch_to_str_only_near_a_placed_member() {
        // Two placed members on R-A (GD 0 to each other) + one STR-only GD 2 away → suggested R-A.
        // A distant member (GD 12) forms its own cluster, not suggested R-A.
        let members = vec![
            m("kit1", Some("R-A"), &haplo(0)),
            m("kit2", Some("R-A"), &haplo(0)),
            m("kit3", None, &haplo(2)),  // STR-only, close → R-A
            m("kit4", None, &haplo(12)), // STR-only, far → own cluster, no suggestion
        ];
        let c = cluster_ystr(&members, &ClusterOpts::default());

        // The R-A cluster holds kit1/kit2 (placed) + kit3 (suggested).
        let ra = c
            .clusters
            .iter()
            .find(|cl| cl.branch.as_deref() == Some("R-A"))
            .unwrap();
        assert_eq!(ra.confirmed_count(), 2);
        assert_eq!(ra.suggested_count(), 1);
        let k3 = ra.members.iter().find(|x| x.label == "kit3").unwrap();
        let sug = k3.suggested.as_ref().unwrap();
        assert_eq!(sug.branch, "R-A");
        assert_eq!(sug.gd, 2);
        assert!(sug.confidence > 0.5);

        // kit4 is in a different cluster (GD 12 from the others, above the 8/100 link threshold).
        let other = c
            .clusters
            .iter()
            .find(|cl| cl.members.iter().any(|x| x.label == "kit4"))
            .unwrap();
        assert!(
            other.members.iter().all(|x| x.label != "kit1"),
            "kit4 must not join the R-A cluster"
        );
        // kit4 has no placed cluster-mate → no suggestion.
        let k4 = other.members.iter().find(|x| x.label == "kit4").unwrap();
        assert!(k4.suggested.is_none());
    }

    #[test]
    fn too_few_markers_are_unclustered() {
        let members = vec![
            m("placed", Some("R-A"), &haplo(0)),
            m("tiny", None, &[("DYS393", "13"), ("DYS390", "24")]), // 2 markers < min
        ];
        let c = cluster_ystr(&members, &ClusterOpts::default());
        assert_eq!(c.unclustered.len(), 1);
        assert_eq!(c.unclustered[0].label, "tiny");
        assert!(c.unclustered[0].suggested.is_none());
    }
}
