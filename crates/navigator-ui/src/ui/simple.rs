//! Simple mode's subject view: a left section rail and the dedicated panels it switches between.
//!
//! Simple mode began as a single vertical scroll of every brief section stacked in a row. That grew
//! past the point where a casual reader could find anything in it — the paternal line, the ancestry
//! donut, the relatives list and the chat all lived in one column several screens tall. This module
//! splits that column into panels reached from a persistent rail, so the vertical extent of any one
//! screen is bounded and the rail itself doubles as the summary (each item carries its own headline
//! value: the terminal haplogroup, the top ancestry, the match count).
//!
//! The panels are ordered as a story rather than by subsystem: [`SimplePanel::Story`] is the landing
//! synopsis, the two lineage panels reach the furthest back in time, [`SimplePanel::Ancestry`] runs
//! deep origins → recent populations in chronological order, and [`SimplePanel::Relatives`] lands in
//! the present day. Nothing here computes anything: it renders the precomputed [`SubjectBrief`] plus
//! the live (online) relative suggestions, exactly as the old scroll did.

use super::*;

/// Width of the section rail. Wide enough for a two-line item (label over its value) at the default
/// text size without wrapping the haplogroup names, which are the longest values it carries.
const RAIL_WIDTH: f32 = 208.0;

/// Fixed size of one at-a-glance tile. Fixed rather than content-sized so `horizontal_wrapped` can
/// reflow the row (see [`NavigatorApp::simple_glance_grid`]); tall enough for a two-line value.
const TILE_SIZE: egui::Vec2 = egui::vec2(186.0, 76.0);

impl NavigatorApp {
    /// The whole Simple-mode subject body: the rail on the left, the selected panel on the right.
    ///
    /// The reference-download and "not analyzed yet" prompts render *above* the split rather than in
    /// a panel — both block every panel equally, so burying either one behind a rail click would let
    /// a user wander an empty view without being told why it's empty.
    pub(crate) fn simple_subject_view(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.separator();
        self.reference_prompt(ui);
        self.simple_analyze_prompt(ui, guid);
        ui.add_space(4.0);

        let full_h = ui.available_height();
        ui.horizontal_top(|ui| {
            ui.allocate_ui_with_layout(
                egui::vec2(RAIL_WIDTH, full_h),
                egui::Layout::top_down(egui::Align::Min),
                |ui| self.simple_rail(ui, guid),
            );
            ui.separator();
            ui.vertical(|ui| {
                egui::ScrollArea::vertical()
                    .id_salt("simple_panel_scroll")
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        ui.add_space(2.0);
                        match self.simple_panel {
                            SimplePanel::Story => self.simple_story_panel(ui, guid),
                            SimplePanel::Paternal => self.simple_lineage_panel(ui, guid, LineageKind::Paternal),
                            SimplePanel::Maternal => self.simple_lineage_panel(ui, guid, LineageKind::Maternal),
                            SimplePanel::Ancestry => self.simple_ancestry_panel(ui, guid),
                            SimplePanel::Relatives => self.simple_relatives_panel(ui, guid),
                            SimplePanel::Test => self.simple_test_panel(ui, guid),
                        }
                        ui.add_space(16.0);
                    });
            });
        });
    }

    // ---------------------------------------------------------------------------------------------
    // The rail
    // ---------------------------------------------------------------------------------------------

    /// The section rail: one clickable item per panel, each showing its own headline value so the
    /// rail reads as a summary of the whole subject without opening anything.
    fn simple_rail(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        // Values first (immutable reads), so the item loop can mutate `simple_panel` freely.
        let items: Vec<(SimplePanel, &'static str, &'static str, Option<String>)> = SimplePanel::ALL
            .iter()
            .map(|(panel, icon, key)| (*panel, *icon, self.tr(key), self.simple_rail_value(*panel, guid)))
            .collect();
        let unavailable = self.tr("simple.rail.noData");

        let mut pick = None;
        for (panel, icon, label, value) in &items {
            let selected = self.simple_panel == *panel;
            let fill = if selected {
                ui.visuals().selection.bg_fill.gamma_multiply(0.45)
            } else {
                egui::Color32::TRANSPARENT
            };
            let resp = egui::Frame::none()
                .fill(fill)
                .rounding(egui::Rounding::same(6.0))
                .inner_margin(egui::Margin::symmetric(10.0, 8.0))
                .show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new(*icon).weak());
                        ui.label(if selected {
                            egui::RichText::new(*label).strong()
                        } else {
                            egui::RichText::new(*label)
                        });
                    });
                    // A missing value is stated, not hidden: an item with no line under it reads as
                    // a rendering gap, where "Not available yet" reads as a fact about the data.
                    if let Some(v) = value {
                        ui.label(egui::RichText::new(v).weak().small());
                    } else if *panel != SimplePanel::Story {
                        ui.label(egui::RichText::new(unavailable).weak().small().italics());
                    }
                })
                .response
                .interact(egui::Sense::click())
                .on_hover_cursor(egui::CursorIcon::PointingHand);
            if resp.clicked() {
                pick = Some(*panel);
            }
            ui.add_space(2.0);
        }
        if let Some(panel) = pick {
            self.simple_panel = panel;
        }

        // The bridge to the full power-user view lives at the foot of the rail, where it's reachable
        // from every panel rather than only from the bottom of one long scroll.
        ui.add_space(12.0);
        ui.separator();
        ui.add_space(6.0);
        if ui
            .link(self.tr("brief.seeData"))
            .on_hover_text(self.tr("brief.seeDataHint"))
            .clicked()
        {
            self.enter_advanced_mode();
        }
    }

    /// The one-line value the rail shows under a panel's name, or `None` when that panel has nothing
    /// yet (the rail then says so explicitly).
    ///
    /// Gated on the brief belonging to `guid`: the brief is rebuilt asynchronously after a subject
    /// switch, and an ungated read would spend those frames labelling the new person with the
    /// previous one's haplogroups.
    fn simple_rail_value(&self, panel: SimplePanel, guid: SampleGuid) -> Option<String> {
        let brief = self.subject_brief.as_ref().filter(|(g, _)| *g == guid).map(|(_, b)| b);
        match panel {
            SimplePanel::Story => brief.map(|b| b.headline.test_chip.clone()),
            SimplePanel::Paternal => brief.and_then(|b| b.paternal.as_ref()).map(|l| l.haplogroup.clone()),
            SimplePanel::Maternal => brief.and_then(|b| b.maternal.as_ref()).map(|l| l.haplogroup.clone()),
            SimplePanel::Ancestry => brief
                .and_then(|b| b.ancestry.as_ref())
                .and_then(|a| a.super_populations.first())
                .map(|sp| format!("{:.0}% {}", sp.percentage, sp.super_population)),
            SimplePanel::Relatives => {
                let n = self.simple_relative_count();
                (n > 0).then(|| format!("{n} {}", self.tr("simple.relatives.count")))
            }
            SimplePanel::Test => brief.map(|b| b.test.quality_phrase.clone()),
        }
    }

    /// How many distinct people the Relatives panel would list: network suggestions plus any
    /// completed exchange that has no suggestion behind it (a confirmed relative is still a relative
    /// after the suggestion that introduced them has aged out of the AppView's list).
    fn simple_relative_count(&self) -> usize {
        let suggested: std::collections::HashSet<&str> = self
            .ibd_suggestions
            .iter()
            .map(|s| s.suggested_sample_guid.as_str())
            .collect();
        let extra = self
            .exchange_results
            .iter()
            .filter(|e| !e.partner_sample_ref.as_deref().is_some_and(|r| suggested.contains(r)))
            .count();
        suggested.len() + extra
    }

    // ---------------------------------------------------------------------------------------------
    // Panel: Your story (landing)
    // ---------------------------------------------------------------------------------------------

    /// The landing synopsis: who this person is in one or two sentences, an at-a-glance grid that
    /// jumps into the other panels, and the results chat.
    ///
    /// When the local AI assistant has narrated the brief, that narration *replaces* the one-line
    /// synopsis as the story — it says the same thing at more length and in better prose. The
    /// structured one-liner stays reachable under "Plain summary" so the model's version is never
    /// the only account of the data on the screen.
    fn simple_story_panel(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let Some(brief) = self.subject_brief.as_ref().filter(|(g, _)| *g == guid).map(|(_, b)| b) else {
            self.simple_brief_placeholder(ui);
            return;
        };
        let (name, test_chip, summary) = (
            brief.headline.name.clone(),
            brief.headline.test_chip.clone(),
            brief.headline.summary.clone(),
        );

        ui.heading(&name);
        ui.add_space(2.0);
        chip(
            ui,
            &test_chip,
            egui::Color32::from_rgb(40, 52, 70),
            egui::Color32::from_rgb(150, 190, 240),
        );
        ui.add_space(12.0);

        self.simple_synopsis_card(ui, guid, &summary);
        ui.add_space(12.0);
        self.simple_glance_grid(ui, guid);
        ui.add_space(12.0);
        self.simple_chat_section(ui, guid);
    }

    /// The synopsis card — the AI narration when there is one, the structured one-liner otherwise,
    /// plus the generate/regenerate control when the local assistant is enabled.
    fn simple_synopsis_card(&mut self, ui: &mut egui::Ui, guid: SampleGuid, summary: &str) {
        // Prefer the live stream while generating; fall back to the finalized narration.
        let live = self
            .narration_stream
            .as_ref()
            .filter(|(g, _)| *g == guid)
            .map(|(_, t)| (t.clone(), None));
        let finalized = self
            .brief_narration
            .as_ref()
            .filter(|(g, _)| *g == guid)
            .map(|(_, n)| (n.prose.clone(), Some(n.model.clone())));
        let story = live.or(finalized).filter(|(p, _)| !p.trim().is_empty());

        let title = if story.is_some() {
            self.tr("brief.aiStory")
        } else {
            self.tr("simple.story.title")
        };
        let (disclaimer, model_label, plain_label) = (
            self.tr("brief.aiDisclaimer"),
            self.tr("brief.aiModel"),
            self.tr("simple.story.plainSummary"),
        );
        let ai_enabled = self.ai_enabled;
        let narrating = self.narrating;
        let regen = self.tr("brief.aiRegenerate");
        let polish = self.tr("brief.polishAi");
        let working = self.tr("brief.aiWorking");

        let mut start_narration = false;
        card(ui, title, |ui| {
            match &story {
                Some((prose, model)) => {
                    for para in prose.split("\n\n").filter(|p| !p.trim().is_empty()) {
                        ui.label(para.trim());
                        ui.add_space(6.0);
                    }
                    ui.separator();
                    ui.label(egui::RichText::new(disclaimer).weak().small());
                    if let Some(model) = model {
                        ui.label(egui::RichText::new(format!("{model_label} {model}")).weak().small());
                    }
                    ui.add_space(4.0);
                    egui::CollapsingHeader::new(plain_label)
                        .id_salt("simple_plain_summary")
                        .show(ui, |ui| {
                            ui.label(summary);
                        });
                }
                None => {
                    ui.label(summary);
                }
            }
            if ai_enabled {
                ui.add_space(8.0);
                ui.horizontal(|ui| {
                    let label = if story.is_some() { regen } else { polish };
                    if ui.add_enabled(!narrating, egui::Button::new(label)).clicked() {
                        start_narration = true;
                    }
                    if narrating {
                        ui.spinner();
                        ui.label(egui::RichText::new(working).weak().small());
                    }
                });
            }
        });

        if start_narration {
            self.narrating = true;
            self.brief_narration = None;
            self.narration_stream = Some((guid, String::new())); // live buffer
            let _ = self.tx.send(Command::NarrateBrief(guid));
        }
    }

    /// The at-a-glance grid: one tile per remaining panel, each showing its headline value and
    /// opening that panel when clicked. This is the landing screen's index — it is why the story
    /// panel doesn't need to restate the paternal line, the ancestry donut, and the match list.
    fn simple_glance_grid(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let tiles: Vec<(SimplePanel, &'static str, &'static str, Option<String>)> = SimplePanel::ALL
            .iter()
            .filter(|(p, _, _)| *p != SimplePanel::Story)
            .map(|(panel, icon, key)| (*panel, *icon, self.tr(key), self.simple_rail_value(*panel, guid)))
            .collect();
        let empty = self.tr("simple.rail.noData");
        let heading = self.tr("simple.glance");

        ui.label(egui::RichText::new(heading).strong().size(15.0));
        ui.add_space(6.0);
        let mut pick = None;
        // Each tile is allocated at a fixed size *before* its frame draws. `Frame::show` reserves no
        // space up front, so a frame built straight into `horizontal_wrapped` gives the layout no
        // size to test against and the row never wraps — it just runs off the right edge. Allocating
        // first is what makes the reflow work.
        ui.horizontal_wrapped(|ui| {
            for (panel, icon, label, value) in &tiles {
                let resp = ui
                    .allocate_ui(TILE_SIZE, |ui| {
                        egui::Frame::group(ui.style())
                            .fill(ui.visuals().faint_bg_color)
                            .rounding(egui::Rounding::same(8.0))
                            .inner_margin(egui::Margin::symmetric(14.0, 10.0))
                            .show(ui, |ui| {
                                ui.set_min_size(ui.available_size());
                                ui.vertical(|ui| {
                                    ui.label(egui::RichText::new(format!("{icon}  {label}")).weak().small());
                                    ui.add_space(4.0);
                                    // Values run from "U5a1b1g" to "high-quality (28× average
                                    // depth)", so they wrap inside the tile rather than widening it.
                                    match value {
                                        Some(v) => ui.label(egui::RichText::new(v).size(15.0).strong()),
                                        None => ui.label(egui::RichText::new(empty).size(14.0).weak().italics()),
                                    };
                                });
                            });
                    })
                    .response
                    .interact(egui::Sense::click())
                    .on_hover_cursor(egui::CursorIcon::PointingHand);
                if resp.clicked() {
                    pick = Some(*panel);
                }
            }
        });
        if let Some(panel) = pick {
            self.simple_panel = panel;
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Panel: paternal / maternal line
    // ---------------------------------------------------------------------------------------------

    /// One lineage panel. Renders the brief's lineage card (haplogroup, age, origin, story,
    /// confidence, descent trail) or, when this subject has no such line, an explanation of *why*
    /// there isn't one — "no data" with no reason is the complaint this redesign exists to fix.
    fn simple_lineage_panel(&mut self, ui: &mut egui::Ui, guid: SampleGuid, kind: LineageKind) {
        let Some(brief) = self.subject_brief.as_ref().filter(|(g, _)| *g == guid).map(|(_, b)| b) else {
            self.simple_brief_placeholder(ui);
            return;
        };
        let (lineage, title, dna) = match kind {
            LineageKind::Paternal => (brief.paternal.clone(), self.tr("brief.paternalLine"), DnaType::Y),
            LineageKind::Maternal => (brief.maternal.clone(), self.tr("brief.maternalLine"), DnaType::Mt),
        };
        let Some(lineage) = lineage else {
            let (t, h) = match kind {
                LineageKind::Paternal => ("simple.noPaternal", "simple.noPaternalHint"),
                LineageKind::Maternal => ("simple.noMaternal", "simple.noMaternalHint"),
            };
            empty_state(ui, self.tr(t), self.tr(h));
            return;
        };
        // The descent trail inside the card reads the cached variant profile; ask for it first.
        self.ensure_descent(guid, dna);
        card(ui, title, |ui| self.brief_lineage_card(ui, guid, &lineage));
    }

    // ---------------------------------------------------------------------------------------------
    // Panel: ancestry (deep origins → recent populations)
    // ---------------------------------------------------------------------------------------------

    /// The autosomal ancestry panel, ordered chronologically rather than by method: the ancient
    /// source populations first, then the archaic (Neanderthal) trace, then the continental and
    /// fine-grained populations of the last few thousand years, then the parent-split painting, then
    /// the runs-of-homozygosity read on the two parental lines meeting again.
    ///
    /// Reading it top to bottom is meant to be reading forwards in time. Sections with no data are
    /// skipped silently — each depends on a different optional analysis, and an absent one is a step
    /// not yet run, not a finding.
    fn simple_ancestry_panel(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let Some(brief) = self.subject_brief.as_ref().filter(|(g, _)| *g == guid).map(|(_, b)| b) else {
            self.simple_brief_placeholder(ui);
            return;
        };
        let (ancestry, roh, archaic) = (brief.ancestry.clone(), brief.roh.clone(), brief.archaic.clone());
        let Some(a) = ancestry else {
            empty_state(ui, self.tr("simple.noAncestry"), self.tr("simple.noAncestryHint"));
            return;
        };

        // --- Deep origins: the ancient source populations, furthest back -------------------------
        if !a.ancient_pops.is_empty() {
            simple_era_divider(ui, self.tr("simple.era.deep"));
            let intro = self.tr("brief.ancientIntro");
            let gloss = self.tr("glossary.ancient");
            card(ui, self.tr("brief.ancient"), |ui| {
                ui.label(egui::RichText::new(intro).weak().small()).on_hover_text(gloss);
                ui.add_space(8.0);
                let slices: Vec<(f64, egui::Color32)> = a
                    .ancient_pops
                    .iter()
                    .map(|c| (c.percentage, parse_hex_color(&c.color)))
                    .collect();
                let top_pct = a.ancient_pops.first().map(|c| c.percentage);
                ui.horizontal(|ui| {
                    draw_color_donut(ui, &slices, top_pct);
                    ui.add_space(8.0);
                    ui.vertical(|ui| {
                        for c in &a.ancient_pops {
                            ui.horizontal(|ui| {
                                let (rect, _) = ui.allocate_exact_size(egui::vec2(11.0, 11.0), egui::Sense::hover());
                                ui.painter().rect_filled(rect, 2.0, parse_hex_color(&c.color));
                                ui.label(format!("{}  —  {:.1}%", c.name, c.percentage));
                            });
                        }
                    });
                });
                for c in &a.ancient_pops {
                    if let Some(blurb) = &c.blurb {
                        ui.add_space(6.0);
                        ui.label(egui::RichText::new(&c.name).strong().small());
                        ui.label(egui::RichText::new(blurb).small());
                    }
                }
            });
            ui.add_space(10.0);
        }

        // --- Archaic: the oldest trace of all, from before modern humans left Africa --------------
        if let Some(arch) = &archaic {
            let gloss = self.tr("glossary.archaic");
            card(ui, self.tr("brief.archaic"), |ui| {
                ui.heading(&arch.pattern).on_hover_text(gloss);
                ui.add_space(4.0);
                ui.label(&arch.summary_phrase);
                ui.add_space(4.0);
                // The count, framed as copies-of-copies-assayed exactly as the Advanced card does —
                // never a "percent Neanderthal" (archaic-ancestry design S1/S7).
                ui.label(format!(
                    "{} of {} marker copies",
                    arch.total_copies, arch.possible_copies
                ));
                if let (Some(p), Some(c)) = (arch.percentile, &arch.cohort) {
                    ui.label(
                        egui::RichText::new(format!("More than {p:.0}% of {c} reference samples"))
                            .weak()
                            .small(),
                    );
                }
                ui.add_space(6.0);
                self.ai_explain(ui, guid, SignalKind::Archaic);
            });
            ui.add_space(10.0);
        }

        // --- Recent: the continental and fine-grained populations --------------------------------
        simple_era_divider(ui, self.tr("simple.era.recent"));
        let ancestry_gloss = self.tr("glossary.ancestry");
        let detail_title = self.tr("simple.ancestry.detailTitle");
        card(ui, self.tr("simple.ancestry.whereTitle"), |ui| {
            ui.heading(&a.summary_phrase).on_hover_text(ancestry_gloss);
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                draw_ancestry_donut(ui, &a.super_populations);
                ui.add_space(8.0);
                ui.vertical(|ui| {
                    for sp in &a.super_populations {
                        if sp.percentage >= 0.5 {
                            ui.label(format!("{}  —  {:.1}%", sp.super_population, sp.percentage));
                        }
                    }
                });
            });
            if let Some(interp) = &a.interpretation {
                ui.add_space(8.0);
                ui.label(interp);
            }
            ui.add_space(6.0);
            ui.label(egui::RichText::new(&a.method_note).weak().small());
        });
        ui.add_space(10.0);

        // The fine breakdown gets its own card rather than a collapsed row inside the one above:
        // the panel has the vertical room the old single scroll didn't, and it is the section most
        // readers came for.
        if !a.fine_pops.is_empty() {
            card(ui, detail_title, |ui| {
                for (name, pct) in a.fine_pops.iter().filter(|(_, p)| *p >= 0.5) {
                    ui.horizontal(|ui| {
                        ui.label(name);
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            ui.label(egui::RichText::new(format!("{pct:.1}%")).strong());
                        });
                    });
                }
            });
            ui.add_space(10.0);
        }

        // --- The two parental sides, and where they meet -----------------------------------------
        self.simple_dna_sides_section(ui, guid);
        ui.add_space(10.0);

        if let Some(r) = &roh {
            let roh_gloss = self.tr("glossary.roh");
            card(ui, self.tr("brief.roh"), |ui| {
                ui.heading(&r.pattern).on_hover_text(roh_gloss);
                ui.add_space(4.0);
                ui.label(&r.summary_phrase);
                ui.add_space(4.0);
                ui.label(egui::RichText::new(format!("F_ROH {:.4}", r.f_roh)).weak().small());
                ui.add_space(6.0);
                self.ai_explain(ui, guid, SignalKind::Roh);
            });
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Panel: relatives
    // ---------------------------------------------------------------------------------------------

    /// Genetic relatives, grouped by how close the relationship looks. Two kinds of row land here:
    /// a **confirmed** relative, where an encrypted segment exchange has actually run and produced a
    /// shared-cM total, and a **candidate**, where the AppView has only scored the signals we both
    /// published. They are grouped by the same three bands, but the band is chosen from measured cM
    /// where we have it and from the suggestion's tier where we don't — and the row says which,
    /// because "close family" inferred from a score is a much weaker claim than the same words
    /// backed by 1,800 shared centimorgans.
    fn simple_relatives_panel(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        ui.heading(self.tr("brief.relatives"));
        ui.add_space(4.0);
        ui.label(egui::RichText::new(self.tr("brief.relativesNote")).weak().small());
        ui.add_space(10.0);

        if self.account.is_none() {
            // The heading is already on screen, so this is the hint alone rather than a full
            // empty-state block repeating it.
            ui.label(self.tr("brief.relativesSignIn"));
            return;
        }

        let rows = self.simple_relative_rows();
        let loading = self.loading_ibd_suggestions;
        let mut do_find = false;
        ui.horizontal(|ui| {
            if ui
                .add_enabled(!loading, egui::Button::new(self.tr("brief.relativesFind")))
                .clicked()
            {
                do_find = true;
            }
            if loading {
                ui.spinner();
            }
            if !rows.is_empty() {
                ui.label(
                    egui::RichText::new(format!("{} {}", rows.len(), self.tr("simple.relatives.count")))
                        .weak()
                        .small(),
                );
            }
        });
        ui.add_space(10.0);

        if rows.is_empty() {
            if !loading {
                ui.label(egui::RichText::new(self.tr("brief.relativesEmpty")).weak());
            }
        } else {
            let mut introduce = None;
            for (tier, title_key, note_key) in RelativeTier::ALL {
                let group: Vec<&RelativeRow> = rows.iter().filter(|r| r.tier == tier).collect();
                if group.is_empty() {
                    continue;
                }
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new(self.tr(title_key)).strong().size(14.0));
                    ui.label(egui::RichText::new(format!("({})", group.len())).weak().small());
                });
                ui.label(egui::RichText::new(self.tr(note_key)).weak().small());
                ui.add_space(6.0);
                for row in group {
                    if let Some(g) = self.simple_relative_row(ui, row) {
                        introduce = Some(g);
                    }
                }
                ui.add_space(12.0);
            }
            // Resolve back to the full suggestion: the introduction records both AppView sample
            // handles in the matching ledger, and only the suggestion carries them.
            if let Some(suggestion) = introduce.and_then(|g| {
                self.ibd_suggestions
                    .iter()
                    .find(|s| s.suggested_sample_guid == g)
                    .cloned()
            }) {
                self.status = self.tr("network.introducing").to_string();
                let _ = self.tx.send(Command::RequestIntroduction {
                    suggestion,
                    biosample_guid: Some(guid),
                });
            }
            ui.add_space(4.0);
            self.ai_explain(ui, guid, SignalKind::Ibd);
        }

        if do_find {
            self.loading_ibd_suggestions = true;
            self.status = self.tr("network.finding").to_string();
            let _ = self.tx.send(Command::LoadIbdSuggestions);
        }
    }

    /// One relative row. Returns the sample guid when its Connect button was pressed.
    fn simple_relative_row(&self, ui: &mut egui::Ui, row: &RelativeRow) -> Option<String> {
        let mut introduce = None;
        egui::Frame::group(ui.style())
            .fill(ui.visuals().faint_bg_color)
            .rounding(egui::Rounding::same(6.0))
            .inner_margin(egui::Margin::symmetric(12.0, 10.0))
            .show(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new(row.dots).monospace())
                        .on_hover_text(self.tr(row.strength_key));
                    ui.label(egui::RichText::new(&row.handle).monospace())
                        .on_hover_text(&row.full_id);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        match (&row.intro, &row.sample_guid) {
                            (Some(status), _) => {
                                ui.label(egui::RichText::new(status).weak().small());
                            }
                            (None, Some(g)) => {
                                if ui.button(self.tr("brief.relativesConnect")).clicked() {
                                    introduce = Some(g.clone());
                                }
                            }
                            (None, None) => {}
                        }
                    });
                });
                // The evidence line. Measured sharing beats a score, so it's what's shown when we
                // have it; the tier note above already said which kind of row this is.
                match &row.shared {
                    Some(s) => {
                        ui.label(
                            egui::RichText::new(format!(
                                "{:.0} {} · {} {}",
                                s.total_cm,
                                self.tr("simple.relatives.cmUnit"),
                                s.segments,
                                self.tr("simple.relatives.segUnit")
                            ))
                            .small(),
                        );
                        ui.label(
                            egui::RichText::new(format!(
                                "{} {}",
                                self.tr("simple.relatives.confirmed"),
                                s.relationship
                            ))
                            .weak()
                            .small(),
                        );
                    }
                    None => {
                        if !row.signals.is_empty() {
                            ui.label(
                                egui::RichText::new(format!(
                                    "{} {}",
                                    self.tr("brief.relativesWhy"),
                                    row.signals.join(", ")
                                ))
                                .small(),
                            );
                        }
                        ui.label(
                            egui::RichText::new(self.tr("simple.relatives.estimate"))
                                .weak()
                                .small()
                                .italics(),
                        );
                    }
                }
            });
        ui.add_space(4.0);
        introduce
    }

    /// Build the relative rows: every network suggestion, joined to a completed exchange where one
    /// exists, plus any completed exchange with no suggestion behind it.
    fn simple_relative_rows(&self) -> Vec<RelativeRow> {
        // Completed exchanges, keyed by the partner sample they attested to.
        let by_ref: std::collections::HashMap<&str, &navigator_app::StoredIbdExchange> = self
            .exchange_results
            .iter()
            .filter_map(|e| e.partner_sample_ref.as_deref().map(|r| (r, e)))
            .collect();

        let mut rows: Vec<RelativeRow> = self
            .ibd_suggestions
            .iter()
            .map(|s| {
                let shared = by_ref
                    .get(s.suggested_sample_guid.as_str())
                    .map(|e| SharedDna::from(*e));
                RelativeRow::new(
                    &s.suggested_sample_guid,
                    Some(s.suggested_sample_guid.clone()),
                    s.strength(),
                    s.signals.clone(),
                    self.ibd_intros.get(&s.suggested_sample_guid).cloned(),
                    shared,
                    s.score,
                )
            })
            .collect();

        // Confirmed relatives whose suggestion is gone (or never existed): still relatives.
        let suggested: std::collections::HashSet<&str> = self
            .ibd_suggestions
            .iter()
            .map(|s| s.suggested_sample_guid.as_str())
            .collect();
        for e in &self.exchange_results {
            let id = e.partner_sample_ref.as_deref().unwrap_or(e.partner_did.as_str());
            if e.partner_sample_ref.as_deref().is_some_and(|r| suggested.contains(r)) {
                continue;
            }
            rows.push(RelativeRow::new(
                id,
                None,
                MatchStrength::Strong,
                Vec::new(),
                None,
                Some(SharedDna::from(e)),
                1.0,
            ));
        }

        // Strongest first inside each tier: measured sharing, then the composite score.
        rows.sort_by(|a, b| {
            let key = |r: &RelativeRow| (r.shared.as_ref().map(|s| s.total_cm).unwrap_or(0.0), r.score);
            key(b).partial_cmp(&key(a)).unwrap_or(std::cmp::Ordering::Equal)
        });
        rows
    }

    // ---------------------------------------------------------------------------------------------
    // Panel: your test
    // ---------------------------------------------------------------------------------------------

    /// The test the whole brief rests on: what it is, what it can and can't tell, its quality, the
    /// global caveats, how fresh the narrative content is, and the export.
    fn simple_test_panel(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let Some(brief) = self.subject_brief.as_ref().filter(|(g, _)| *g == guid).map(|(_, b)| b) else {
            self.simple_brief_placeholder(ui);
            return;
        };
        let test = brief.test.clone();
        let caveats = brief.caveats.clone();
        let (pack_status, pack_version, enriched) = (brief.pack_status, brief.pack_version.clone(), brief.enriched);

        card(ui, self.tr("brief.yourTest"), |ui| {
            ui.strong(&test.test_name);
            ui.add_space(2.0);
            ui.label(&test.what_it_tells);
            if let Some(lim) = &test.limitations {
                ui.add_space(2.0);
                ui.label(egui::RichText::new(lim).weak().small());
            }
            ui.add_space(8.0);
            let (bg, fg, mark) = if test.quality_ok {
                (
                    egui::Color32::from_rgb(28, 56, 36),
                    egui::Color32::from_rgb(120, 200, 140),
                    "✔",
                )
            } else {
                (
                    egui::Color32::from_rgb(64, 50, 24),
                    egui::Color32::from_rgb(230, 180, 90),
                    "⚠",
                )
            };
            chip(ui, &format!("{mark}  {}", test.quality_phrase), bg, fg);
        });

        if !caveats.is_empty() {
            ui.add_space(10.0);
            card(ui, self.tr("simple.test.caveats"), |ui| {
                for c in &caveats {
                    ui.label(egui::RichText::new(format!("• {c}")).weak().small());
                }
            });
        }

        self.simple_realign_card(ui, guid);

        ui.add_space(10.0);
        self.export_row(ui, &[navigator_app::ExportRequest::SubjectBriefHtml(guid)]);

        // Pack-status footer (how fresh the descriptions are).
        ui.add_space(10.0);
        let pack_key = match pack_status {
            PackStatus::Downloaded => "brief.packLive",
            PackStatus::Cached => "brief.packCached",
            PackStatus::Bundled => "brief.packOffline",
            PackStatus::Unavailable => "brief.packUnavailable",
        };
        let mut footer = self.tr(pack_key).to_string();
        if let Some(v) = &pack_version {
            footer = format!("{footer} · {v}");
        }
        if enriched {
            footer = format!("{footer} · {}", self.tr("brief.enriched"));
        }
        ui.label(egui::RichText::new(footer).weak().small());
    }

    /// The offer to rebuild this person's genome against the complete reference, and the progress of
    /// that job once it is running.
    ///
    /// Lives under "Your test" because that is where this belongs honestly: it is a statement about
    /// the limits of the data, alongside the depth and quality phrases. A reader here should never
    /// need to know what a reference build is — only that the reference itself was unfinished, which
    /// is a fact about the reference and not about their test.
    ///
    /// The copy is careful on two points that are easy to get wrong. The completed assembly is
    /// **genome-wide** — a full autosomal sequence as well as the first complete Y — and Y discovery
    /// is merely what Navigator puts it to work on today; writing it as a paternal-line feature
    /// would misdescribe what T2T actually delivered. And the finished genome is one donor's: of
    /// European ancestry, carrying a J1a Y. Whose DNA a reference represents is exactly the sort of
    /// thing that goes unsaid and should not.
    ///
    /// Whether to offer it at all was decided in the app layer, on the brief; see
    /// `navigator_domain::brief::RealignOffer`.
    fn simple_realign_card(&mut self, ui: &mut egui::Ui, guid: SampleGuid) {
        let offer = self.subject_brief.as_ref().and_then(|(_, b)| b.realign_offer.clone());

        // A run in progress wins over the offer, so the card a user just started reports itself
        // rather than continuing to invite the thing it is already doing.
        //
        // Matched on the *subject*, not the alignment. Comparing alignment ids meant this card had
        // nothing to compare against once the offer was gone — and it went on to claim any running
        // job, so a page open on one person announced another person's realignment as their own.
        let running = self.realign.clone().filter(|state| state.biosample_guid == Some(guid));

        if let Some(state) = running {
            ui.add_space(10.0);
            card(ui, self.tr("simple.realign.title"), |ui| match &state.finished {
                None => {
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label(format!(
                            "{} — {} {} / {}",
                            state.label,
                            self.tr("simple.realign.step"),
                            state.step,
                            state.total
                        ));
                    });
                    if !state.detail.is_empty() {
                        ui.label(egui::RichText::new(&state.detail).weak().small());
                    }
                    ui.add_space(6.0);
                    ui.label(
                        egui::RichText::new(self.tr("simple.realign.runningHint"))
                            .weak()
                            .small(),
                    );
                }
                Some(RealignFinished::Done { .. }) => {
                    ui.label(self.tr("simple.realign.done"));
                }
                Some(RealignFinished::Cancelled) => {
                    ui.label(egui::RichText::new(self.tr("simple.realign.cancelled")).weak());
                }
                Some(RealignFinished::Failed(message)) => {
                    ui.colored_label(ui.visuals().error_fg_color, self.tr("simple.realign.failed"));
                    ui.label(egui::RichText::new(message).weak().small());
                }
            });
            return;
        }

        let Some(offer) = offer else { return };

        ui.add_space(10.0);
        card(ui, self.tr("simple.realign.title"), |ui| {
            ui.label(self.tr("simple.realign.body"));
            ui.add_space(6.0);
            ui.label(self.tr("simple.realign.use"));
            ui.add_space(6.0);
            // Whose genome the finished reference actually is, and whose Y — stated rather than
            // left for a reader to assume it represents everyone equally.
            ui.label(egui::RichText::new(self.tr("simple.realign.note")).weak().small());
            ui.add_space(6.0);
            ui.label(egui::RichText::new(self.tr("simple.realign.cost")).weak().small());
            ui.add_space(8.0);
            // Opens the confirmation rather than starting: see `simple_realign_confirm_modal`.
            if ui.button(self.tr("simple.realign.action")).clicked() {
                self.simple_realign_confirm = Some(offer.clone());
            }
        });
    }

    // ---------------------------------------------------------------------------------------------

    /// Shown by every panel while the brief for this subject is still being built (or when there is
    /// none to build), so no panel renders as blank.
    fn simple_brief_placeholder(&self, ui: &mut egui::Ui) {
        if self.subject_brief_loading {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label(egui::RichText::new(self.tr("brief.building")).weak());
            });
        } else {
            ui.label(egui::RichText::new(self.tr("brief.empty")).weak());
        }
    }
}

/// A muted "you are here in time" divider between the ancestry panel's eras.
fn simple_era_divider(ui: &mut egui::Ui, label: &str) {
    ui.add_space(2.0);
    ui.label(egui::RichText::new(label).weak().small());
    ui.add_space(6.0);
}

/// How close a relative looks. Note this is a *presentation* band, not a relationship estimate: it
/// decides which heading a person is listed under, and each heading states what the band was read
/// from (see [`NavigatorApp::simple_relatives_panel`]).
#[derive(Clone, Copy, PartialEq, Eq)]
enum RelativeTier {
    Close,
    Extended,
    Distant,
}

impl RelativeTier {
    /// `(tier, heading key, what-this-band-means key)` in display order.
    const ALL: [(RelativeTier, &'static str, &'static str); 3] = [
        (
            RelativeTier::Close,
            "simple.relatives.close",
            "simple.relatives.closeNote",
        ),
        (
            RelativeTier::Extended,
            "simple.relatives.extended",
            "simple.relatives.extendedNote",
        ),
        (
            RelativeTier::Distant,
            "simple.relatives.distant",
            "simple.relatives.distantNote",
        ),
    ];

    /// Measured sharing decides the band when there is any; the composite score's tier decides it
    /// otherwise. The cM cutoffs are the conservative ends of the standard bands — 1,200 cM is about
    /// the floor for an aunt/uncle-or-closer relationship, 200 cM about the floor for a range that
    /// still resolves to a nameable cousin rather than "somewhere back there".
    fn of(shared: Option<&SharedDna>, strength: MatchStrength) -> Self {
        match shared {
            Some(s) if s.total_cm >= 1200.0 => RelativeTier::Close,
            Some(s) if s.total_cm >= 200.0 => RelativeTier::Extended,
            Some(_) => RelativeTier::Distant,
            None => match strength {
                MatchStrength::Strong => RelativeTier::Close,
                MatchStrength::Likely => RelativeTier::Extended,
                MatchStrength::Possible => RelativeTier::Distant,
            },
        }
    }
}

/// The measured half of a relative row: what a completed segment exchange actually found.
struct SharedDna {
    total_cm: f64,
    segments: i64,
    relationship: String,
}

impl From<&navigator_app::StoredIbdExchange> for SharedDna {
    fn from(e: &navigator_app::StoredIbdExchange) -> Self {
        SharedDna {
            total_cm: e.total_shared_cm,
            segments: e.segment_count,
            relationship: e.relationship.clone(),
        }
    }
}

/// One row of the relatives panel, pre-resolved so rendering does no interpretation.
struct RelativeRow {
    /// Short pseudonymous handle shown in the list.
    handle: String,
    /// The full identifier, on hover.
    full_id: String,
    /// The sample guid to request an introduction to, when this row came from a suggestion.
    sample_guid: Option<String>,
    tier: RelativeTier,
    dots: &'static str,
    strength_key: &'static str,
    signals: Vec<String>,
    intro: Option<String>,
    shared: Option<SharedDna>,
    score: f64,
}

impl RelativeRow {
    fn new(
        id: &str,
        sample_guid: Option<String>,
        strength: MatchStrength,
        signals: Vec<String>,
        intro: Option<String>,
        shared: Option<SharedDna>,
        score: f64,
    ) -> Self {
        let (dots, strength_key) = match strength {
            MatchStrength::Strong => ("⚫⚫⚫", "brief.matchStrong"),
            MatchStrength::Likely => ("⚫⚫⚪", "brief.matchLikely"),
            MatchStrength::Possible => ("⚫⚪⚪", "brief.matchPossible"),
        };
        RelativeRow {
            handle: id.chars().take(12).collect(),
            full_id: id.to_string(),
            sample_guid,
            tier: RelativeTier::of(shared.as_ref(), strength),
            dots,
            strength_key,
            signals,
            intro,
            shared,
            score,
        }
    }
}
