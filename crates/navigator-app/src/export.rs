//! Result exports (gap §6): pure formatters that turn a cached analysis result into a shareable
//! file body — TSV / HTML / BED. Kept free of I/O and `App` so they're trivially unit-testable; the
//! app layer loads the result and writes the returned `String` to the user-chosen path.

use navigator_analysis::coverage::CoverageResult;
use navigator_analysis::haplo::CallState;
use navigator_analysis::ibd::IbdSegment;
use navigator_analysis::mtvariants::MtVariant;
use navigator_analysis::read_metrics::ReadMetrics;
use navigator_domain::ancestry::AncestryResult;
use navigator_domain::brief::{LineageBrief, SubjectBrief};
use navigator_domain::reconciliation::DnaType;

use crate::{BranchReport, DescentReport};

/// Minimal HTML text escaping for the small, controlled strings we embed (population names etc.).
fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}

/// Shared inline stylesheet for the HTML exports (self-contained — no external assets).
const HTML_STYLE: &str = "body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:2rem;color:#222}\
h1{font-size:1.3rem}h2{font-size:1rem;margin-top:1.5rem}\
table{border-collapse:collapse;margin-top:.5rem}\
th,td{border:1px solid #ddd;padding:4px 10px;text-align:right;font-variant-numeric:tabular-nums}\
th{background:#f3f3f3}td:first-child,th:first-child{text-align:left}\
.bar{height:12px;background:#3b82f6;border-radius:2px;display:inline-block;vertical-align:middle}\
.meta{color:#666;font-size:.85rem}";

// ---- coverage ----------------------------------------------------------------

/// Coverage as TSV: a `#`-commented genome-wide metrics header, then a per-contig table (joining the
/// samtools-style depth stats with the GATK-style callable breakdown).
pub fn coverage_tsv(cov: &CoverageResult) -> String {
    let mut out = String::new();
    out.push_str("# DUNavigator coverage export\n");
    out.push_str(&format!("# genome_territory\t{}\n", cov.genome_territory));
    out.push_str(&format!("# mean_coverage\t{:.4}\n", cov.mean_coverage));
    out.push_str(&format!("# median_coverage\t{:.1}\n", cov.median_coverage));
    out.push_str(&format!("# sd_coverage\t{:.4}\n", cov.sd_coverage));
    out.push_str(&format!("# mad_coverage\t{:.4}\n", cov.mad_coverage));
    out.push_str(&format!("# callable_bases\t{}\n", cov.callable_bases));
    out.push_str(&format!("# pct_exc_mapq\t{:.4}\n", cov.pct_exc_mapq));
    out.push_str(&format!("# pct_exc_baseq\t{:.4}\n", cov.pct_exc_baseq));
    out.push_str(&format!(
        "# pct_1x\t{:.2}\tpct_5x\t{:.2}\tpct_10x\t{:.2}\tpct_15x\t{:.2}\tpct_20x\t{:.2}\tpct_30x\t{:.2}\tpct_50x\t{:.2}\n",
        cov.pct_1x, cov.pct_5x, cov.pct_10x, cov.pct_15x, cov.pct_20x, cov.pct_30x, cov.pct_50x
    ));
    out.push_str(
        "contig\tlength\treads\tcovered_bases\tpct_covered\tmean_depth\tmean_baseq\tmean_mapq\t\
         callable\tno_coverage\tlow_coverage\texcessive\tpoor_mapq\tref_n\n",
    );
    for s in &cov.contig_coverage_stats {
        let c = cov.contig_callable.iter().find(|m| m.contig == s.contig);
        let (callable, nocov, low, exc, poor, refn) = c.map_or((0, 0, 0, 0, 0, 0), |m| {
            (
                m.callable,
                m.no_coverage,
                m.low_coverage,
                m.excessive_coverage,
                m.poor_mapping_quality,
                m.ref_n,
            )
        });
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            s.contig,
            s.end_pos,
            s.num_reads,
            s.cov_bases,
            s.coverage,
            s.mean_depth,
            s.mean_base_q,
            s.mean_map_q,
            callable,
            nocov,
            low,
            exc,
            poor,
            refn,
        ));
    }
    out
}

/// Coverage as a self-contained HTML page (genome-wide summary + per-contig table).
pub fn coverage_html(cov: &CoverageResult, label: &str) -> String {
    let mut rows = String::new();
    for s in &cov.contig_coverage_stats {
        let c = cov.contig_callable.iter().find(|m| m.contig == s.contig);
        let callable = c.map_or(0, |m| m.callable);
        rows.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.2}</td><td>{:.2}</td><td>{}</td></tr>",
            esc(&s.contig),
            s.end_pos,
            s.num_reads,
            s.coverage,
            s.mean_depth,
            callable,
        ));
    }
    format!(
        "<!doctype html><html><head><meta charset=\"utf-8\"><title>Coverage — {title}</title>\
         <style>{style}</style></head><body>\
         <h1>Coverage — {title}</h1>\
         <table>\
         <tr><th>Metric</th><th>Value</th></tr>\
         <tr><td>Genome territory</td><td>{territory}</td></tr>\
         <tr><td>Mean coverage</td><td>{mean:.2}×</td></tr>\
         <tr><td>Median coverage</td><td>{median:.0}×</td></tr>\
         <tr><td>SD coverage</td><td>{sd:.2}</td></tr>\
         <tr><td>MAD coverage</td><td>{mad:.2}</td></tr>\
         <tr><td>Callable bases</td><td>{callable}</td></tr>\
         <tr><td>% ≥1×</td><td>{p1:.2}</td></tr>\
         <tr><td>% ≥10×</td><td>{p10:.2}</td></tr>\
         <tr><td>% ≥20×</td><td>{p20:.2}</td></tr>\
         <tr><td>% ≥30×</td><td>{p30:.2}</td></tr>\
         </table>\
         <h2>Per-contig</h2>\
         <table><tr><th>Contig</th><th>Length</th><th>Reads</th><th>% covered</th><th>Mean depth</th><th>Callable</th></tr>{rows}</table>\
         </body></html>",
        title = esc(label),
        style = HTML_STYLE,
        territory = cov.genome_territory,
        mean = cov.mean_coverage,
        median = cov.median_coverage,
        sd = cov.sd_coverage,
        mad = cov.mad_coverage,
        callable = cov.callable_bases,
        p1 = cov.pct_1x,
        p10 = cov.pct_10x,
        p20 = cov.pct_20x,
        p30 = cov.pct_30x,
    )
}

// ---- read metrics ------------------------------------------------------------

/// Read metrics as a two-column `metric<TAB>value` TSV (Picard-style metric names).
pub fn read_metrics_tsv(m: &ReadMetrics) -> String {
    let mut out = String::from("# DUNavigator read-metrics export\nmetric\tvalue\n");
    let mut row = |k: &str, v: String| out.push_str(&format!("{k}\t{v}\n"));
    row("TOTAL_READS", m.total_reads.to_string());
    row("PF_READS", m.pf_reads.to_string());
    row("PF_READS_ALIGNED", m.pf_reads_aligned.to_string());
    row("PCT_PF_READS_ALIGNED", format!("{:.4}", m.pct_pf_reads_aligned));
    row("READS_ALIGNED_IN_PAIRS", m.reads_aligned_in_pairs.to_string());
    row(
        "PCT_READS_ALIGNED_IN_PAIRS",
        format!("{:.4}", m.pct_reads_aligned_in_pairs),
    );
    row("PROPER_PAIRS", m.proper_pairs.to_string());
    row("PCT_PROPER_PAIRS", format!("{:.4}", m.pct_proper_pairs));
    row("MEDIAN_READ_LENGTH", format!("{:.1}", m.median_read_length));
    row("MEAN_READ_LENGTH", format!("{:.2}", m.mean_read_length));
    row("SD_READ_LENGTH", format!("{:.2}", m.std_read_length));
    row("MIN_READ_LENGTH", m.min_read_length.to_string());
    row("MAX_READ_LENGTH", m.max_read_length.to_string());
    row("MEDIAN_INSERT_SIZE", format!("{:.1}", m.median_insert_size));
    row("MEAN_INSERT_SIZE", format!("{:.2}", m.mean_insert_size));
    row("SD_INSERT_SIZE", format!("{:.2}", m.std_insert_size));
    row("MIN_INSERT_SIZE", m.min_insert_size.to_string());
    row("MAX_INSERT_SIZE", m.max_insert_size.to_string());
    row("PAIR_ORIENTATION", m.pair_orientation.as_str().to_string());
    row("PCT_CHIMERAS", format!("{:.4}", m.pct_chimeras));
    row("MEAN_MAPPING_QUALITY", format!("{:.2}", m.mean_mapping_quality));
    out
}

// ---- ancestry ----------------------------------------------------------------

/// Ancestry as a single TSV table — super-population then fine-population rows distinguished by a
/// `level` column — under a `#`-commented metadata header.
pub fn ancestry_tsv(a: &AncestryResult) -> String {
    let mut out = String::new();
    out.push_str("# DUNavigator ancestry export\n");
    out.push_str(&format!("# method\t{}\n", a.method));
    out.push_str(&format!("# panel_type\t{}\n", a.panel_type));
    out.push_str(&format!(
        "# snps_analyzed\t{}\tsnps_with_genotype\t{}\tsnps_missing\t{}\n",
        a.snps_analyzed, a.snps_with_genotype, a.snps_missing
    ));
    out.push_str(&format!("# confidence_level\t{:.4}\n", a.confidence_level));
    if let Some(fd) = a.fit_distance {
        out.push_str(&format!("# fit_distance\t{fd:.4}\n"));
    }
    out.push_str(&format!(
        "# pipeline_version\t{}\treference_version\t{}\n",
        a.pipeline_version, a.reference_version
    ));
    out.push_str("level\tcode\tname\tpercentage\trank\tci_lower\tci_upper\n");
    for s in &a.super_population_summary {
        out.push_str(&format!(
            "super\t{}\t{}\t{:.2}\t\t\t\n",
            s.super_population, s.super_population, s.percentage
        ));
    }
    for c in &a.components {
        out.push_str(&format!(
            "population\t{}\t{}\t{:.2}\t{}\t{:.2}\t{:.2}\n",
            c.population_code,
            c.population_name,
            c.percentage,
            c.rank,
            c.confidence_interval.lower,
            c.confidence_interval.upper,
        ));
    }
    out
}

/// Ancestry as a self-contained HTML page: metadata, super-population bars, and the fine-population
/// table with confidence intervals.
pub fn ancestry_html(a: &AncestryResult) -> String {
    let mut bars = String::new();
    for s in &a.super_population_summary {
        bars.push_str(&format!(
            "<div><span class=\"bar\" style=\"width:{w}px\"></span> {name} — {pct:.1}%</div>",
            w = (s.percentage * 2.0).round().max(1.0) as i64,
            name = esc(&s.super_population),
            pct = s.percentage,
        ));
    }
    let mut rows = String::new();
    for c in &a.components {
        rows.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.2}</td><td>{:.1} – {:.1}</td></tr>",
            c.rank,
            esc(&c.population_code),
            esc(&c.population_name),
            c.percentage,
            c.confidence_interval.lower,
            c.confidence_interval.upper,
        ));
    }
    let call_rate = if a.snps_analyzed > 0 {
        100.0 * a.snps_with_genotype as f64 / a.snps_analyzed as f64
    } else {
        0.0
    };
    format!(
        "<!doctype html><html><head><meta charset=\"utf-8\"><title>Ancestry — {method}</title>\
         <style>{style}</style></head><body>\
         <h1>Ancestry estimate</h1>\
         <p class=\"meta\">Method {method} · panel {panel} · {wg}/{an} SNPs genotyped ({cr:.1}% call rate) · confidence {conf:.0}%</p>\
         <h2>Continental ancestry</h2>{bars}\
         <h2>Populations</h2>\
         <table><tr><th>Rank</th><th>Code</th><th>Population</th><th>%</th><th>95% CI</th></tr>{rows}</table>\
         </body></html>",
        method = esc(&a.method),
        panel = esc(&a.panel_type),
        style = HTML_STYLE,
        wg = a.snps_with_genotype,
        an = a.snps_analyzed,
        cr = call_rate,
        conf = a.confidence_level * 100.0,
    )
}

// ---- mtDNA variants ----------------------------------------------------------

/// mtDNA variants (vs rCRS) as TSV: position, compact notation, region, ref/alt, type.
/// The Y-DNA / mtDNA **descent report** (root→terminal lineage) as TSV: one row per defining SNP of
/// each node on the path, with the subject's call state and observed base. Mirrors the on-screen
/// per-node descent grid so it can be shared / diffed outside the app.
/// TSV for a [`BranchReport`] — one row per defining marker in the reported subtree, with the
/// sample's observed base + call state + evidence. Shareable for placement spot-checks / researcher
/// exchange. Missing evidence renders as `.` (VCF convention).
pub fn branch_report_tsv(report: &BranchReport) -> String {
    let dna = match report.dna {
        DnaType::Y => "Y-DNA",
        DnaType::Mt => "mtDNA",
    };
    let state = |s: CallState| match s {
        CallState::Derived => "derived",
        CallState::Ancestral => "ancestral",
        CallState::NoCall => "nocall",
    };
    let gt = |s: CallState| match s {
        CallState::Derived => "1",
        CallState::Ancestral => "0",
        CallState::NoCall => ".",
    };
    let u = |o: Option<u32>| o.map(|v| v.to_string()).unwrap_or_else(|| ".".to_string());
    let ad = |o: Option<(u32, u32)>| o.map(|(r, a)| format!("{r},{a}")).unwrap_or_else(|| ".".to_string());
    let base = |o: Option<char>| o.map(|c| c.to_string()).unwrap_or_else(|| ".".to_string());

    let (d, a, n) = report.counts();
    let mut out = format!(
        "# DUNavigator {dna} branch report — node {} ({}); {d} derived / {a} ancestral / {n} no-call\n",
        report.root, report.contig
    );
    out.push_str("node\tparent\tmarker\tchrom\tpos\tancestral\tderived\tobserved_base\tstatus\tGT\tAD\tDP\tGQ\tsource\tnote\n");
    for r in &report.rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.node,
            r.parent,
            r.marker,
            report.contig,
            r.position,
            r.ancestral,
            r.derived,
            base(r.observed_base),
            state(r.state),
            gt(r.state),
            ad(r.ad),
            u(r.dp),
            u(r.gq),
            r.source,
            r.note,
        ));
    }
    out
}

pub fn descent_tsv(report: &DescentReport) -> String {
    let dna = match report.dna {
        DnaType::Y => "Y-DNA",
        DnaType::Mt => "mtDNA",
    };
    let state = |s: CallState| match s {
        CallState::Derived => "derived",
        CallState::Ancestral => "ancestral",
        CallState::NoCall => "nocall",
    };
    let mut out = format!("# DUNavigator {dna} descent report — terminal {}\n", report.terminal);
    out.push_str("haplogroup\tterminal\tsnp\tposition\tancestral\tderived\tstate\tobserved_base\n");
    for node in &report.nodes {
        for snp in &node.snps {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                node.name,
                if node.is_terminal { "yes" } else { "" },
                snp.name,
                snp.position,
                snp.ancestral,
                snp.derived,
                state(snp.state),
                snp.base.map(|c| c.to_string()).unwrap_or_default(),
            ));
        }
    }
    out
}

pub fn mtdna_variants_tsv(variants: &[MtVariant]) -> String {
    let mut out = String::from("# DUNavigator mtDNA variants vs rCRS (NC_012920.1)\n");
    out.push_str("position\tnotation\tregion\tref\talt\ttype\n");
    for v in variants {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{:?}\n",
            v.position,
            v.notation(),
            v.region().label(),
            v.reference,
            v.alternate,
            v.kind,
        ));
    }
    out
}

// ---- IBD segments ------------------------------------------------------------

/// IBD segments as TSV (`chromosome  start  end  length_cm  snp_count`), 1-based bp. The match
/// browser's "Export segments CSV" — a tab-delimited table for downstream analysis / sharing.
pub fn ibd_segments_tsv(segments: &[IbdSegment]) -> String {
    let mut out = String::from("# DUNavigator IBD segments export\n");
    out.push_str("chromosome\tstart_position\tend_position\tlength_cm\tsnp_count\thalf_identical\n");
    for s in segments {
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.2}\t{}\t{}\n",
            s.chromosome,
            s.start_position,
            s.end_position,
            s.length_cm,
            s.snp_count.map(|n| n.to_string()).unwrap_or_default(),
            s.is_half_identical.map(|b| b.to_string()).unwrap_or_default(),
        ));
    }
    out
}

// ---- callable loci -----------------------------------------------------------

/// Callable intervals as BED4 (`contig  start  end  CALLABLE`), 0-based half-open. `per_contig` is
/// `(contig, [(start0, end)])` as produced by `coverage::callable_intervals`.
pub fn callable_bed(per_contig: &[(String, Vec<(i64, i64)>)]) -> String {
    let mut out = String::from("# DUNavigator callable loci (CALLABLE state) — BED4, 0-based half-open\n");
    for (contig, intervals) in per_contig {
        for (start, end) in intervals {
            out.push_str(&format!("{contig}\t{start}\t{end}\tCALLABLE\n"));
        }
    }
    out
}

// ---- subject brief ("DNA Story") ---------------------------------------------

/// One lineage block of the brief HTML.
fn lineage_html(lb: &LineageBrief, title: &str) -> String {
    let mut s = format!("<h2>{}</h2>\n<p class=\"hg\">{}</p>\n", esc(title), esc(&lb.haplogroup));
    if let Some(anc) = &lb.matched_ancestor {
        s.push_str(&format!(
            "<p class=\"meta\">The description below is for {}, an ancestor on this line.</p>\n",
            esc(anc)
        ));
    }
    if let Some(age) = &lb.age_phrase {
        s.push_str(&format!("<p><em>{}</em></p>\n", esc(age)));
    }
    if let Some(origin) = &lb.origin_phrase {
        s.push_str(&format!("<p>{}</p>\n", esc(origin)));
    }
    if let Some(story) = &lb.story {
        s.push_str(&format!("<p>{}</p>\n", esc(story)));
    }
    s.push_str(&format!("<p class=\"meta\">{}</p>\n", esc(&lb.confidence_phrase)));
    if !lb.sources.is_empty() {
        s.push_str(&format!(
            "<p class=\"meta\">Source: {}</p>\n",
            esc(&lb.sources.join(", "))
        ));
    }
    s
}

/// The subject brief as a self-contained "DNA Story" HTML document — the casual-reader report a user
/// can save or print. Mirrors the Simple-mode card stack. When an AI narration is provided (a cached
/// "Polish with AI" result), it leads the document as a clearly-labelled, additive section above the
/// structured facts.
pub fn subject_brief_html(b: &SubjectBrief, narration: Option<&crate::NarratedBrief>) -> String {
    let mut body = String::new();
    body.push_str(&format!("<h1>{} — Your DNA Story</h1>\n", esc(&b.headline.name)));
    body.push_str(&format!("<p class=\"meta\">{}</p>\n", esc(&b.headline.test_chip)));
    body.push_str(&format!("<p>{}</p>\n", esc(&b.headline.summary)));

    if let Some(n) = narration {
        body.push_str("<h2>Your DNA Story (AI-assisted)</h2>\n");
        for para in n.prose.split("\n\n").filter(|p| !p.trim().is_empty()) {
            body.push_str(&format!("<p>{}</p>\n", esc(para.trim())));
        }
        body.push_str(&format!(
            "<p class=\"meta\">AI-assisted from your results (model: {}) — the verified facts follow below.</p>\n",
            esc(&n.model)
        ));
    }

    if let Some(p) = &b.paternal {
        body.push_str(&lineage_html(p, "Your paternal line (Y-DNA)"));
    }
    if let Some(m) = &b.maternal {
        body.push_str(&lineage_html(m, "Your maternal line (mtDNA)"));
    }

    if let Some(a) = &b.ancestry {
        body.push_str("<h2>Your ancestry</h2>\n");
        body.push_str(&format!("<p class=\"hg\">{}</p>\n", esc(&a.summary_phrase)));
        body.push_str("<table><tr><th>Population</th><th>Share</th></tr>\n");
        for sp in a.super_populations.iter().filter(|s| s.percentage >= 0.5) {
            body.push_str(&format!(
                "<tr><td>{}</td><td>{:.1}%</td></tr>\n",
                esc(&sp.super_population),
                sp.percentage
            ));
        }
        body.push_str("</table>\n");
        if let Some(interp) = &a.interpretation {
            body.push_str(&format!("<p>{}</p>\n", esc(interp)));
        }
        if !a.fine_pops.is_empty() {
            body.push_str("<h3>Detailed populations</h3>\n<table><tr><th>Population</th><th>Share</th></tr>\n");
            for (name, pct) in a.fine_pops.iter().filter(|(_, p)| *p >= 0.5) {
                body.push_str(&format!("<tr><td>{}</td><td>{:.1}%</td></tr>\n", esc(name), pct));
            }
            body.push_str("</table>\n");
        }
        if !a.ancient_pops.is_empty() {
            body.push_str("<h3>Ancient ancestry</h3>\n");
            body.push_str("<table><tr><th>Component</th><th>Share</th></tr>\n");
            for c in &a.ancient_pops {
                body.push_str(&format!(
                    "<tr><td>{}</td><td>{:.1}%</td></tr>\n",
                    esc(&c.name),
                    c.percentage
                ));
            }
            body.push_str("</table>\n");
            for c in &a.ancient_pops {
                if let Some(blurb) = &c.blurb {
                    body.push_str(&format!("<p><strong>{}.</strong> {}</p>\n", esc(&c.name), esc(blurb)));
                }
            }
        }
        body.push_str(&format!("<p class=\"meta\">{}</p>\n", esc(&a.method_note)));
    }

    if let Some(r) = &b.roh {
        body.push_str("<h2>Shared ancestry</h2>\n");
        body.push_str(&format!("<p class=\"hg\">{}</p>\n", esc(&r.pattern)));
        body.push_str(&format!("<p>{}</p>\n", esc(&r.summary_phrase)));
        body.push_str(&format!("<p class=\"meta\">F_ROH {:.4}</p>\n", r.f_roh));
    }

    body.push_str("<h2>Your test</h2>\n");
    body.push_str(&format!("<p class=\"hg\">{}</p>\n", esc(&b.test.test_name)));
    body.push_str(&format!("<p>{}</p>\n", esc(&b.test.what_it_tells)));
    if let Some(lim) = &b.test.limitations {
        body.push_str(&format!("<p class=\"meta\">{}</p>\n", esc(lim)));
    }
    let mark = if b.test.quality_ok { "✓" } else { "⚠" };
    body.push_str(&format!("<p>{mark} {}</p>\n", esc(&b.test.quality_phrase)));

    if !b.caveats.is_empty() {
        body.push_str("<h2>Notes</h2>\n<ul>\n");
        for c in &b.caveats {
            body.push_str(&format!("<li>{}</li>\n", esc(c)));
        }
        body.push_str("</ul>\n");
    }

    let mut footer = match b.pack_status {
        navigator_domain::brief::PackStatus::Downloaded => "Descriptions updated online".to_string(),
        navigator_domain::brief::PackStatus::Cached => "Descriptions from your last online update".to_string(),
        navigator_domain::brief::PackStatus::Bundled => "Offline descriptions".to_string(),
        navigator_domain::brief::PackStatus::Unavailable => "Descriptions unavailable".to_string(),
    };
    if let Some(v) = &b.pack_version {
        footer.push_str(&format!(" · {v}"));
    }
    if b.enriched {
        footer.push_str(" · live haplogroup data");
    }
    body.push_str(&format!("<p class=\"meta\">{}</p>\n", esc(&footer)));

    format!(
        "<!doctype html><html><head><meta charset=\"utf-8\"><title>{} — DNA Story</title>\
<style>{HTML_STYLE}\nh3{{font-size:.95rem;margin-top:1rem}}.hg{{font-size:1.1rem;font-weight:600;margin:.2rem 0}}\
td,th{{text-align:left}}</style></head><body>\n{body}</body></html>",
        esc(&b.headline.name)
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_analysis::mtvariants::{MtVariant, MtVariantKind};
    use navigator_domain::ancestry::{AncestryResult, ConfidenceInterval, PopulationComponent, SuperPopulationSummary};

    #[test]
    fn subject_brief_html_renders_sections() {
        use navigator_domain::brief::{
            AncestryBrief, AncientComponent, Headline, LineageBrief, LineageKind, PackStatus, SubjectBrief, TestBrief,
        };
        let brief = SubjectBrief {
            headline: Headline {
                name: "James".into(),
                test_chip: "Whole Genome Sequencing".into(),
                summary: "Your data places James's paternal line at R-FGC29071.".into(),
            },
            paternal: Some(LineageBrief {
                kind: LineageKind::Paternal,
                haplogroup: "R-FGC29071".into(),
                lineage_path: vec!["R".into(), "R-M269".into(), "R-FGC29071".into()],
                matched_ancestor: Some("R-M269".into()),
                age_phrase: Some("formed roughly 6,400 years ago".into()),
                origin_phrase: Some("associated with the steppe".into()),
                story: Some("A common Western European lineage.".into()),
                confidence_phrase: "strong placement, confirmed across multiple tests".into(),
                sources: vec!["YFull".into()],
            }),
            maternal: None,
            ancestry: Some(AncestryBrief {
                summary_phrase: "Predominantly European".into(),
                super_populations: vec![SuperPopulationSummary {
                    super_population: "European".into(),
                    percentage: 98.0,
                    populations: vec![],
                }],
                fine_pops: vec![("British".into(), 60.0)],
                ancient_pops: vec![AncientComponent {
                    code: "Steppe".into(),
                    name: "Steppe pastoralist".into(),
                    percentage: 50.0,
                    color: "#4e79a7".into(),
                    blurb: Some("Bronze Age steppe migrants.".into()),
                }],
                interpretation: Some("European ancestry spans the continent.".into()),
                method_note: "estimated from 400,000 genome-wide markers".into(),
            }),
            roh: None,
            test: TestBrief {
                test_name: "Whole Genome Sequencing".into(),
                what_it_tells: "Reads your whole genome.".into(),
                limitations: None,
                quality_phrase: "high-quality (30× average depth)".into(),
                quality_ok: true,
            },
            needs_analysis: false,
            caveats: vec![],
            pack_version: Some("2026.06-seed".into()),
            pack_status: PackStatus::Bundled,
            enriched: true,
        };
        let narration = crate::NarratedBrief {
            prose: "You are mostly European.".into(),
            model: "test-model".into(),
        };
        let html = subject_brief_html(&brief, Some(&narration));
        assert!(html.starts_with("<!doctype html>"));
        assert!(html.contains("Your DNA Story"));
        assert!(html.contains("AI-assisted"));
        assert!(html.contains("You are mostly European."));
        assert!(html.contains("R-FGC29071"));
        assert!(html.contains("Steppe pastoralist"));
        assert!(html.contains("Predominantly European"));
        // No narration → no AI section.
        assert!(!subject_brief_html(&brief, None).contains("AI-assisted"));
        assert!(html.contains("live haplogroup data"));
    }

    #[test]
    fn coverage_tsv_has_header_and_joined_contig_rows() {
        use navigator_analysis::coverage::{ContigCallableMetrics, ContigCoverageStats, CoverageResult};
        let cov = CoverageResult {
            genome_territory: 1000,
            mean_coverage: 28.4,
            callable_bases: 900,
            contig_coverage_stats: vec![ContigCoverageStats {
                contig: "chr1".into(),
                start_pos: 1,
                end_pos: 500,
                num_reads: 42,
                cov_bases: 480,
                coverage: 96.0,
                mean_depth: 30.1,
                mean_base_q: 35.0,
                mean_map_q: 58.0,
                histogram: vec![],
            }],
            contig_callable: vec![ContigCallableMetrics {
                contig: "chr1".into(),
                ref_n: 5,
                callable: 470,
                no_coverage: 10,
                low_coverage: 15,
                excessive_coverage: 0,
                poor_mapping_quality: 0,
            }],
            ..Default::default()
        };
        let tsv = coverage_tsv(&cov);
        assert!(tsv.contains("# mean_coverage\t28.4000"));
        assert!(tsv.contains("contig\tlength\treads"));
        // contig row joins stats (num_reads=42) with callable (callable=470).
        assert!(tsv
            .lines()
            .any(|l| l.starts_with("chr1\t500\t42\t480\t96.00\t30.10\t35.0\t58.0\t470\t10\t15\t0\t0\t5")));
        // HTML variant renders without panicking and includes the title.
        assert!(coverage_html(&cov, "KANE-0001").contains("Coverage — KANE-0001"));
    }

    fn ancestry() -> AncestryResult {
        AncestryResult {
            method: "ADMIXTURE".into(),
            panel_type: "aims".into(),
            snps_analyzed: 18537,
            snps_with_genotype: 17900,
            snps_missing: 637,
            components: vec![PopulationComponent {
                population_code: "GBR".into(),
                population_name: "British".into(),
                percentage: 87.4,
                confidence_interval: ConfidenceInterval {
                    lower: 85.0,
                    upper: 89.8,
                },
                rank: 1,
            }],
            super_population_summary: vec![SuperPopulationSummary {
                super_population: "EUR".into(),
                percentage: 100.0,
                populations: vec!["GBR".into()],
            }],
            confidence_level: 0.96,
            fit_distance: None,
            pipeline_version: "1".into(),
            reference_version: "chm13v2".into(),
            pca_coordinates: None,
        }
    }

    #[test]
    fn ancestry_tsv_distinguishes_super_and_population_rows() {
        let tsv = ancestry_tsv(&ancestry());
        assert!(tsv.contains("# method\tADMIXTURE"));
        assert!(tsv.lines().any(|l| l.starts_with("super\tEUR\tEUR\t100.00")));
        assert!(tsv
            .lines()
            .any(|l| l.starts_with("population\tGBR\tBritish\t87.40\t1\t85.00\t89.80")));
        // HTML escapes + renders.
        let html = ancestry_html(&ancestry());
        assert!(html.contains("Method ADMIXTURE"));
        assert!(html.contains("British"));
    }

    #[test]
    fn mtdna_tsv_uses_notation_and_region() {
        let vs = vec![
            MtVariant {
                position: 263,
                reference: "A".into(),
                alternate: "G".into(),
                kind: MtVariantKind::Substitution,
            },
            MtVariant {
                position: 16519,
                reference: "T".into(),
                alternate: "C".into(),
                kind: MtVariantKind::Substitution,
            },
        ];
        let tsv = mtdna_variants_tsv(&vs);
        assert!(tsv
            .lines()
            .any(|l| l.starts_with("263\t263A>G\tHVR2\tA\tG\tSubstitution")));
        assert!(tsv
            .lines()
            .any(|l| l.starts_with("16519\t16519T>C\tHVR1\tT\tC\tSubstitution")));
    }

    #[test]
    fn callable_bed_is_bed4_half_open() {
        let bed = callable_bed(&[("chr1".into(), vec![(0, 1000), (2000, 2500)])]);
        assert!(bed.lines().any(|l| l == "chr1\t0\t1000\tCALLABLE"));
        assert!(bed.lines().any(|l| l == "chr1\t2000\t2500\tCALLABLE"));
    }

    #[test]
    fn ibd_segments_tsv_has_header_and_rows() {
        let tsv = ibd_segments_tsv(&[IbdSegment {
            chromosome: "chr1".into(),
            start_position: 1,
            end_position: 10_000_000,
            length_cm: 12.5,
            snp_count: Some(80),
            is_half_identical: Some(false),
        }]);
        assert!(tsv.lines().next().unwrap().starts_with('#'));
        assert!(tsv.contains("chromosome\tstart_position\tend_position\tlength_cm\tsnp_count\thalf_identical"));
        assert!(tsv.lines().any(|l| l == "chr1\t1\t10000000\t12.50\t80\tfalse"));
    }

    #[test]
    fn descent_tsv_has_header_and_per_snp_rows() {
        use navigator_analysis::haplo::{NodeEvidence, SnpEvidence};
        let report = DescentReport {
            dna: DnaType::Y,
            terminal: "R-FGC29071".into(),
            nodes: vec![NodeEvidence {
                name: "R-M269".into(),
                is_terminal: false,
                snps: vec![
                    SnpEvidence {
                        name: "M269".into(),
                        position: 22739367,
                        ancestral: "T".into(),
                        derived: "C".into(),
                        state: CallState::Derived,
                        base: Some('C'),
                    },
                    SnpEvidence {
                        name: "L23".into(),
                        position: 6753511,
                        ancestral: "G".into(),
                        derived: "A".into(),
                        state: CallState::NoCall,
                        base: None,
                    },
                ],
            }],
        };
        let tsv = descent_tsv(&report);
        assert!(tsv.lines().next().unwrap().starts_with("# DUNavigator Y-DNA"));
        assert!(tsv.contains("haplogroup\tterminal\tsnp\tposition\tancestral\tderived\tstate\tobserved_base"));
        assert!(tsv.lines().any(|l| l == "R-M269\t\tM269\t22739367\tT\tC\tderived\tC"));
        assert!(tsv.lines().any(|l| l == "R-M269\t\tL23\t6753511\tG\tA\tnocall\t"));
    }

    #[test]
    fn branch_report_tsv_has_header_and_evidence_rows() {
        use crate::BranchRow;
        let report = BranchReport {
            dna: DnaType::Y,
            root: "R-FGC29071".into(),
            contig: "chrY".into(),
            gvcf_backed: true,
            rows: vec![
                BranchRow {
                    node: "R-FGC29071".into(),
                    parent: "R-FGC29067".into(),
                    marker: "FGC29069".into(),
                    position: 14583465,
                    ancestral: "G".into(),
                    derived: "T".into(),
                    observed_base: Some('T'),
                    state: CallState::Derived,
                    ad: Some((0, 11)),
                    dp: Some(11),
                    gq: Some(99),
                    source: "gvcf_variant",
                    note: String::new(),
                },
                BranchRow {
                    node: "R-MF41134".into(),
                    parent: "R-FGC29071".into(),
                    marker: "MF41134".into(),
                    position: 12803849,
                    ancestral: "C".into(),
                    derived: "T".into(),
                    observed_base: Some('C'),
                    state: CallState::Ancestral,
                    ad: None,
                    dp: None,
                    gq: Some(99),
                    source: "gvcf_refblock",
                    note: "hom-ref block".into(),
                },
            ],
        };
        let tsv = branch_report_tsv(&report);
        assert!(tsv.lines().next().unwrap().starts_with("# DUNavigator Y-DNA branch report — node R-FGC29071"));
        assert!(tsv.contains("node\tparent\tmarker\tchrom\tpos\tancestral\tderived\tobserved_base\tstatus\tGT\tAD\tDP\tGQ\tsource\tnote"));
        assert!(tsv
            .lines()
            .any(|l| l == "R-FGC29071\tR-FGC29067\tFGC29069\tchrY\t14583465\tG\tT\tT\tderived\t1\t0,11\t11\t99\tgvcf_variant\t"));
        // Ref-block row: AD/DP omitted (.), GQ kept, ancestral, hom-ref note.
        assert!(tsv
            .lines()
            .any(|l| l == "R-MF41134\tR-FGC29071\tMF41134\tchrY\t12803849\tC\tT\tC\tancestral\t0\t.\t.\t99\tgvcf_refblock\thom-ref block"));
    }
}
