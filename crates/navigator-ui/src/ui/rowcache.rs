//! Per-frame view caches for the three tables whose display rows are expensive to derive.
//!
//! egui is immediate-mode: a render fn runs on **every** frame, up to 60 times a second, for as long
//! as its tab is open. Deriving display rows inline therefore repeats that work at frame rate even
//! when nothing changed — 6 `String` clones and a natural sort per subject, 15 per project member, or
//! (worst) a full scan of an autosomal consensus panel of ~1.2M sites just to count how many match
//! the filter.
//!
//! Each cache here records the inputs its rows were derived from and rebuilds only when one of them
//! differs. The data inputs are nearly all written in `NavigatorApp::drain_events`, so they collapse
//! to a single `data_epoch`; the rest are UI state the user changes directly (selection, language,
//! sort column, filter text). Rebuilding is always safe — the failure mode to avoid is *not*
//! rebuilding — so `data_epoch` is bumped on every event rather than per affected field.
//!
//! Every key also carries the **length** of the collection the rows were derived from. These caches
//! store indices, so a stale index into a shrunken collection would panic rather than merely look
//! wrong — and not every mutation goes through an event (`select_project`, for one, clears
//! `project_report` directly on a click). Keying on the length makes a cached index into a
//! different-length collection impossible to observe, without having to enumerate every mutator.

use super::*;
use crate::i18n::Lang;

/// One row of the subjects table: the subject it selects, plus its rendered cell text.
pub(crate) struct SubjectRow {
    pub(crate) guid: SampleGuid,
    pub(crate) cells: [String; 6],
}

/// Filtered + sorted rows for the subjects table.
#[derive(Default)]
pub(crate) struct SubjectRowCache {
    /// `None` until the first build.
    lang: Option<Lang>,
    epoch: u64,
    /// `all_biosamples.len()` when the rows were built — see the module docs.
    len: usize,
    selected: Option<SampleGuid>,
    sort: Option<usize>,
    ascending: bool,
    filters: Vec<String>,
    pub(crate) rows: Vec<SubjectRow>,
}

/// Row order for the project Report table. The 15 cells of display text per member are what the
/// filter and the natural sort run over, but they are only ever needed during a rebuild — the body
/// renders from the live report rows (the lite badge and action buttons need more than flat text),
/// so `order` indexes `NavigatorApp::project_report` and the text itself is dropped.
#[derive(Default)]
pub(crate) struct ReportRowCache {
    built: bool,
    epoch: u64,
    /// `project_report.len()` when `order` was built — see the module docs.
    len: usize,
    sort: Option<usize>,
    ascending: bool,
    filters: Vec<String>,
    pub(crate) order: Vec<usize>,
}

/// Indices of the variants matching a variant table's status filter and search text.
///
/// The autosomal consensus panel runs to ~1.2M sites and a WGS Y profile to thousands. Only the first
/// `CAP` matches are ever rendered, but the *count* ("N of M matching") needs the whole scan — so
/// without this the entire panel was walked on every frame the tab was open.
#[derive(Default)]
pub(crate) struct VariantRows {
    built: bool,
    epoch: u64,
    /// `variants.len()` when `rows` was built — see the module docs.
    len: usize,
    filter: Option<YVariantStatus>,
    query: String,
    rows: Vec<u32>,
}

impl VariantRows {
    /// Indices into `variants` of the entries `keep` accepts, rescanned only when `epoch`, `filter`,
    /// or `query` changes. `keep` must depend on nothing else that varies between frames.
    pub(crate) fn get<T>(
        &mut self,
        epoch: u64,
        filter: Option<YVariantStatus>,
        query: &str,
        variants: &[T],
        keep: impl Fn(&T) -> bool,
    ) -> &[u32] {
        if !self.built
            || self.epoch != epoch
            || self.len != variants.len()
            || self.filter != filter
            || self.query != query
        {
            self.rows.clear();
            let matching = variants
                .iter()
                .enumerate()
                .filter(|(_, v)| keep(v))
                .map(|(i, _)| i as u32);
            self.rows.extend(matching);
            self.built = true;
            self.epoch = epoch;
            self.len = variants.len();
            self.filter = filter;
            self.query.clear();
            self.query.push_str(query);
        }
        &self.rows
    }
}

impl NavigatorApp {
    /// Rebuild [`Self::subject_rows`] if any input changed. Split from the render fn so the caller can
    /// then borrow `self.subject_rows` and `self.subjects_table_ctl` as disjoint fields.
    pub(crate) fn refresh_subject_rows(&mut self) {
        let c = &self.subject_rows;
        let ctl = &self.subjects_table_ctl;
        let current = c.lang == Some(self.lang)
            && c.epoch == self.data_epoch
            && c.len == self.all_biosamples.len()
            && c.selected == self.selected_sample
            && c.sort == ctl.sort_col()
            && c.ascending == ctl.ascending()
            && c.filters == ctl.filters_raw();
        if current {
            return;
        }

        let mut rows: Vec<SubjectRow> = self
            .all_biosamples
            .iter()
            .map(|s| {
                // Y/mt from the bulk per-subject summary; the selected row prefers the freshly
                // loaded consensus (reflects a just-run assignment before the summary reloads).
                let sel = self.selected_sample == Some(s.guid);
                let summary = self.haplo_summary.get(&s.guid);
                let y = sel
                    .then(|| self.consensus_y.as_ref().map(|c| c.haplogroup.clone()))
                    .flatten()
                    .or_else(|| summary.and_then(|(y, _)| y.clone()))
                    .unwrap_or_else(|| "-".into());
                let mt = sel
                    .then(|| self.consensus_mt.as_ref().map(|c| c.haplogroup.clone()))
                    .flatten()
                    .or_else(|| summary.and_then(|(_, m)| m.clone()))
                    .unwrap_or_else(|| "-".into());
                // Analysis status: Complete once every alignment is analyzed, Pending while any is
                // not (e.g. a just-imported file); a subject with no alignments has no status.
                let status = match self.subject_status.get(&s.guid) {
                    Some(SubjectAnalysisStatus::Complete) => self.tr("subjectStatus.complete"),
                    Some(SubjectAnalysisStatus::Pending) => self.tr("subjectStatus.pending"),
                    None => "-",
                };
                SubjectRow {
                    guid: s.guid,
                    cells: [
                        s.donor_identifier.clone(),
                        y,
                        mt,
                        s.sex.clone().unwrap_or_else(|| "-".into()),
                        s.center_name.clone().unwrap_or_else(|| "-".into()),
                        status.to_string(),
                    ],
                }
            })
            .collect();

        // Inline per-column filters (AND across columns), then natural-sort by the active column.
        for col in 0..SUBJECT_COLS.len() {
            let f = self.subjects_table_ctl.filter_norm(col);
            if !f.is_empty() {
                rows.retain(|r| r.cells[col].to_lowercase().contains(&f));
            }
        }
        if let Some(c) = self.subjects_table_ctl.sort_col() {
            let asc = self.subjects_table_ctl.ascending();
            rows.sort_by(|a, b| {
                let o = natural_cmp(&a.cells[c], &b.cells[c]);
                if asc {
                    o
                } else {
                    o.reverse()
                }
            });
        }

        let ctl = &self.subjects_table_ctl;
        self.subject_rows = SubjectRowCache {
            lang: Some(self.lang),
            epoch: self.data_epoch,
            len: self.all_biosamples.len(),
            selected: self.selected_sample,
            sort: ctl.sort_col(),
            ascending: ctl.ascending(),
            filters: ctl.filters_raw().to_vec(),
            rows,
        };
    }

    /// Rebuild [`Self::report_rows`] if the report or the table controls changed. As with
    /// [`Self::refresh_subject_rows`], kept separate so the caller can borrow the cache and the
    /// controls at once.
    pub(crate) fn refresh_report_rows(&mut self, actions_col: usize) {
        let c = &self.report_rows;
        let ctl = &self.report_table_ctl;
        let current = c.built
            && c.epoch == self.data_epoch
            && c.len == self.project_report.len()
            && c.sort == ctl.sort_col()
            && c.ascending == ctl.ascending()
            && c.filters == ctl.filters_raw();
        if current {
            return;
        }

        // Display text per cell — the basis for inline filtering and natural sort (the body renders
        // from the live report rows so the lite badge + action buttons stay rich).
        let texts: Vec<[String; 15]> = self
            .project_report
            .iter()
            .map(|r| {
                [
                    r.biosample.donor_identifier.clone(),
                    r.alignment_count.to_string(),
                    fmt_depth(r.mean_coverage),
                    fmt_depth(r.median_coverage),
                    fmt_pct(r.pct_10x),
                    fmt_pct(r.pct_20x),
                    r.callable_bases.map(|v| v.to_string()).unwrap_or_else(|| "—".into()),
                    r.y_haplogroup.clone().unwrap_or_else(|| "—".into()),
                    r.mt_haplogroup.clone().unwrap_or_else(|| "—".into()),
                    r.sex.clone().unwrap_or_else(|| "—".into()),
                    fmt_depth(r.mean_read_length),
                    fmt_pct(r.pct_aligned),
                    fmt_depth(r.median_insert_size),
                    r.sv_count.map(|v| v.to_string()).unwrap_or_else(|| "—".into()),
                    String::new(),
                ]
            })
            .collect();

        // Filter (AND across columns) then natural-sort the row order.
        let mut order: Vec<usize> = (0..self.project_report.len()).collect();
        let active_filters: Vec<(usize, String)> = (0..actions_col)
            .filter_map(|c| {
                let f = self.report_table_ctl.filter_norm(c);
                (!f.is_empty()).then_some((c, f))
            })
            .collect();
        if !active_filters.is_empty() {
            order.retain(|&i| {
                active_filters
                    .iter()
                    .all(|(c, f)| texts[i][*c].to_lowercase().contains(f))
            });
        }
        if let Some(c) = self.report_table_ctl.sort_col() {
            if c < actions_col {
                let asc = self.report_table_ctl.ascending();
                order.sort_by(|&a, &b| {
                    let o = natural_cmp(&texts[a][c], &texts[b][c]);
                    if asc {
                        o
                    } else {
                        o.reverse()
                    }
                });
            }
        }

        let ctl = &self.report_table_ctl;
        self.report_rows = ReportRowCache {
            built: true,
            epoch: self.data_epoch,
            len: self.project_report.len(),
            sort: ctl.sort_col(),
            ascending: ctl.ascending(),
            filters: ctl.filters_raw().to_vec(),
            order,
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    /// Count how many entries the predicate was asked about, so "did it rescan?" is observable.
    fn scan(cache: &mut VariantRows, epoch: u64, query: &str, data: &[i32], seen: &Cell<usize>) -> Vec<u32> {
        cache
            .get(epoch, None, query, data, |v| {
                seen.set(seen.get() + 1);
                v % 2 == 0
            })
            .to_vec()
    }

    #[test]
    fn variant_rows_rescan_only_when_an_input_changes() {
        let data: Vec<i32> = (0..6).collect();
        let mut cache = VariantRows::default();
        let seen = Cell::new(0);

        assert_eq!(scan(&mut cache, 1, "", &data, &seen), vec![0, 2, 4]);
        assert_eq!(seen.get(), 6, "first use scans the collection");

        // The whole point: repeated frames with identical inputs must not rescan.
        for _ in 0..30 {
            assert_eq!(scan(&mut cache, 1, "", &data, &seen), vec![0, 2, 4]);
        }
        assert_eq!(seen.get(), 6, "unchanged inputs must not rescan");

        scan(&mut cache, 2, "", &data, &seen);
        assert_eq!(seen.get(), 12, "a new data epoch rescans");
        scan(&mut cache, 2, "rs1", &data, &seen);
        assert_eq!(seen.get(), 18, "new search text rescans");

        let before = seen.get();
        cache.get(2, Some(YVariantStatus::Conflict), "rs1", &data, |v| v % 2 == 0);
        assert!(seen.get() == before, "the status filter is part of the key, not the predicate here");
    }

    /// Indices are only ever handed back for a collection of the length they were derived from —
    /// otherwise a `select_project`-style clear would leave them pointing past the end.
    #[test]
    fn variant_rows_never_index_past_a_shrunken_collection() {
        let data: Vec<i32> = (0..6).collect();
        let mut cache = VariantRows::default();
        let seen = Cell::new(0);
        assert_eq!(scan(&mut cache, 1, "", &data, &seen), vec![0, 2, 4]);

        // Same epoch, same filter, same query — only the collection shrank (mutated outside the
        // event loop). The rows must still be valid indices into it.
        let shrunk: Vec<i32> = vec![4];
        let rows = scan(&mut cache, 1, "", &shrunk, &seen);
        assert!(
            rows.iter().all(|&i| (i as usize) < shrunk.len()),
            "stale indices {rows:?} would panic when used to index {} entries",
            shrunk.len()
        );
        assert_eq!(rows, vec![0]);

        // ...and an emptied collection yields no rows at all.
        assert!(scan(&mut cache, 1, "", &[], &seen).is_empty());
    }
}
