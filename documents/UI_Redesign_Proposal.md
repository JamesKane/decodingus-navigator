# UI Redesign Proposal: DUNavigator

**Date:** December 2024
**Status:** Draft Proposal
**Target Users:** Genetic genealogists and scientists managing large subject collections

---

## Executive Summary

This proposal outlines a modern UI redesign for DUNavigator, transitioning from the current split-pane layout to a **dashboard-centric, entity-focused design**. The redesign prioritizes:

- Efficient management of large subject collections (100+ individuals)
- Rich per-subject analysis views (Y-DNA, mtDNA, Ancestry, IBD)
- Side-by-side comparison capabilities
- Internationalization (i18n) support from the ground up

---

## Current Architecture

### Existing Layout

```
+-----------------------------------------------------+
|            TopBar (Login/Settings)                   |
+----------------------+------------------------------+
|                      |                              |
|   Left Panel         |     Right Panel              |
|   (Navigation)       |     (Details)                |
|                      |                              |
| - Projects List      |  Project / Subject Detail    |
| - Filter             |  - Metadata                  |
| - Add Project Btn    |  - Edit/Delete Actions       |
|                      |  - Haplogroup Info           |
| - Subjects List      |  - Data Tabs                 |
| - Filter             |    - Sequencing              |
| - Add Subject Btn    |    - Chip/Array              |
| - Save Workspace     |    - STR Profiles            |
|                      |                              |
+----------------------+------------------------------+
|  StatusBar (Online/Sync/Conflicts/Cache Info)       |
+-----------------------------------------------------+
```

### Current Pain Points

1. **Flat lists don't scale** - Simple filtered lists become unwieldy at 100+ subjects
2. **Single detail view** - Cannot compare subjects side-by-side
3. **Analysis is buried** - Results hidden in tabs and dialogs
4. **No visual haplogroup navigation** - Text-only tree display
5. **Project/Subject separation** - Hard to see cross-project relationships
6. **No batch operations** - One-at-a-time workflow

---

## Proposed Architecture: Entity-Centric Dashboard

### Design Philosophy

Analysis results and IBD matches are **attributes of subjects**, not independent entities. The navigation should reflect this by having only entity-focused top-level tabs, with analysis and matching living within the subject detail view.

### Top-Level Navigation

```
+-------------------------------------------------------------------------+
|  [=]  Decoding-Us Navigator          [Search...]           [User] [Cfg] |
+-------------------------------------------------------------------------+
|  [Dashboard]  |  [Subjects]  |  [Projects]                              |
+-------------------------------------------------------------------------+
|                                                                         |
|                    [ TAB CONTENT AREA ]                                 |
|                                                                         |
+-------------------------------------------------------------------------+
|  * Online  |  Synced  |  Cache: 2.3GB  |  47 subjects loaded            |
+-------------------------------------------------------------------------+
```

**Three tabs only:**
- **Dashboard** - Aggregate overview, pending work queue
- **Subjects** - The primary workspace (individuals with all their data)
- **Projects** - Research groupings

---

## Tab 1: Dashboard

**Purpose:** At-a-glance overview for quick orientation and pending work management.

```
+-------------------------------------------------------------------------+
|  WORKSPACE OVERVIEW                                                     |
+-------------------------------------------------------------------------+
|                                                                         |
|  WORKSPACE SUMMARY                                                      |
|  +-----------+ +-----------+ +-----------+ +-----------+                |
|  |    47     | |     5     | |    38     | |    12     |                |
|  | Subjects  | | Projects  | | with Y-DNA| |IBD Matches|                |
|  +-----------+ +-----------+ +-----------+ +-----------+                |
|                                                                         |
|  PENDING WORK                                              [Run All]    |
|  +---------------------------------------------------------------------+|
|  | * Jane Doe (DU-002) - Y-DNA analysis pending                        ||
|  | * Bob Johnson (DU-003) - mtDNA analysis pending                     ||
|  | * Alice Wang (DU-004) - Coverage metrics pending                    ||
|  | ! Tom Wilson (DU-005) - Analysis failed (click to retry)            ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  HAPLOGROUP DISTRIBUTIONS                                               |
|  +----------------------------+  +----------------------------+         |
|  |  Y-DNA                     |  |  mtDNA                     |         |
|  |  R1b ============== 17     |  |  H  ================ 14    |         |
|  |  I1  =======        8      |  |  J  =======          6     |         |
|  |  E1b ======         6      |  |  K  ======           5     |         |
|  |  R1a =====          5      |  |  T  =====            4     |         |
|  +----------------------------+  +----------------------------+         |
|                                                                         |
|  RECENT ACTIVITY                                                        |
|  +---------------------------------------------------------------------+|
|  | - 2h ago     John Smith - Haplogroup confirmed R1b-P312             ||
|  | - Yesterday  Jane Doe - Data uploaded                               ||
|  | - Dec 15     Smith Family - 3 members added                         ||
|  +---------------------------------------------------------------------+|
|                                                                         |
+-------------------------------------------------------------------------+
```

---

## Tab 2: Subjects

### Subject Grid View

**Purpose:** Efficiently manage large subject collections with advanced filtering.

```
+-------------------------------------------------------------------------+
| SUBJECTS                                                   [+ Add New]  |
+-------------------------------------------------------------------------+
| [Search name, ID, haplogroup...]     [Filters v]  [Columns v]  [...]    |
+-------------------------------------------------------------------------+
| [ ] | ID      | Name         | Y-DNA      | mtDNA  | Project   | Status |
+-----+---------+--------------+------------+--------+-----------+--------+
| [ ] | DU-001  | John Smith   | R1b-P312   | H1a    | Smith Fam | Done   |
| [ ] | DU-002  | Jane Doe     | -          | J1c    | Doe Res   | Pend   |
| [ ] | DU-003  | Bob Johnson  | I1-M253    | K1a    | Unassigned| Done   |
| [x] | DU-004  | Alice Wang   | O2-M122    | B4     | Asia Proj | Err    |
|     | ...     |              |            |        |           |        |
+-------------------------------------------------------------------------+
| [x] 1 selected                   [Compare] [Batch Analyze] [Add to Proj]|
+-------------------------------------------------------------------------+
```

### Filter Panel (Collapsible)

```
+-------------------------------------------------------------------------+
|                                                                         |
|  Y-DNA Haplogroup        mtDNA Haplogroup        Project                |
|  [R1b v] [Any level v]   [H v] [Any level v]     [All Projects v]       |
|                                                                         |
|  Analysis Status         Data Type               Date Range             |
|  [x] Complete [x] Pend   [x] WGS  [x] Chip       [Last 30 days v]       |
|  [ ] Error   [ ] None    [x] STR  [ ] VCF                               |
|                                                                         |
+-------------------------------------------------------------------------+
```

**Key Features:**
- Multi-select for batch operations
- Column customization (show/hide fields)
- Advanced filtering (haplogroup prefix matching, date ranges)
- Inline status indicators
- Quick compare (select 2+ subjects)

---

## Subject Detail View

When a subject is selected, the right panel (or full-width view) shows a rich detail with **internal tabs** for each analysis domain.

### Header

```
+-------------------------------------------------------------------------+
| <- Back to List                                     [Edit] [Delete] [.].|
+-------------------------------------------------------------------------+
|                                                                         |
|  +--------+   JOHN SMITH                                                |
|  |        |   ID: DU-001 - Male - Smith Family Project                  |
|  |  [Av]  |   Added: Dec 2024                                           |
|  +--------+                                                             |
|                                                                         |
+-------------------------------------------------------------------------+
|  Overview | Y-DNA | mtDNA | Ancestry | IBD Matches | Data Sources       |
+-------------------------------------------------------------------------+
```

### Subject Tab: Overview

Quick summary of all genetic findings:

```
+-------------------------------------------------------------------------+
|  GENETIC SUMMARY                                                        |
|  +---------------------------+  +-----------------------------+         |
|  |  Y-DNA                    |  |  mtDNA                      |         |
|  |  R1b-P312                 |  |  H1a1                       |         |
|  |  Confidence:              |  |  Confidence:                |         |
|  |  ***** HIGH               |  |  ***** HIGH                 |         |
|  |  [View Details ->]        |  |  [View Details ->]          |         |
|  +---------------------------+  +-----------------------------+         |
|                                                                         |
|  +---------------------------+  +-----------------------------+         |
|  |  Ancestry Composition     |  |  IBD Matches                |         |
|  |  o Not yet analyzed       |  |  12 matches found           |         |
|  |                           |  |  Closest: Bob Smith (1847cM)|         |
|  |  [Run Analysis]           |  |  [View All ->]              |         |
|  +---------------------------+  +-----------------------------+         |
|                                                                         |
|  DATA SOURCES                                                           |
|  +---------------------------------------------------------------------+|
|  |  [2 Sequencing Runs]  |  [1 Chip Profile]  |  [1 STR Profile]       ||
|  +---------------------------------------------------------------------+|
+-------------------------------------------------------------------------+
```

### Subject Tab: Y-DNA

All Y-chromosome analysis for this subject:

```
+-------------------------------------------------------------------------+
|  Y-DNA ANALYSIS                                         [Run Analysis]  |
+-------------------------------------------------------------------------+
|                                                                         |
|  TERMINAL HAPLOGROUP                                                    |
|  +---------------------------------------------------------------------+|
|  |                                                                     ||
|  |   R1b-P312                                                          ||
|  |   =====================================================             ||
|  |                                                                     ||
|  |   Phylogenetic Path:                                                ||
|  |   R -> R1 -> R1b -> R1b-M269 -> R1b-L151 -> R1b-P312               ||
|  |                                                                     ||
|  |   Derived: 847    Ancestral: 12                                     ||
|  |   Callable: 98.2%                                                   ||
|  |                                                                     ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  Y-CHROMOSOME IDEOGRAM                                                  |
|  +---------------------------------------------------------------------+|
|  |  [========|====**====|========*=====|=======***====|======]         ||
|  |  p arm                    centromere                    q arm       ||
|  |  * = derived SNP positions                                          ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  SOURCE RECONCILIATION                         Status: Consistent       |
|  +------------------+-------------+-------------+------------------+    |
|  | Source           | Haplogroup  | Derived     | Quality          |    |
|  +------------------+-------------+-------------+------------------+    |
|  | WGS_NovaSeq_2024 | R1b-P312    | 847         | ***** Excellent  |    |
|  | WGS_HiSeq_2022   | R1b-P312    | 612         | ****  Good       |    |
|  | BigY-700 Import  | R1b-P312    | 698         | ****  Good       |    |
|  +------------------+-------------+-------------+------------------+    |
|                                                                         |
|  [View Full Y Profile]  [Export SNP List]  [Compare with Others]        |
|                                                                         |
+-------------------------------------------------------------------------+
```

### Subject Tab: mtDNA

```
+-------------------------------------------------------------------------+
|  mtDNA ANALYSIS                                         [Run Analysis]  |
+-------------------------------------------------------------------------+
|                                                                         |
|  HAPLOGROUP                                                             |
|  +---------------------------------------------------------------------+|
|  |   H1a1                                                              ||
|  |   Path: H -> H1 -> H1a -> H1a1                                      ||
|  |   Confidence: ***** HIGH (156/158 markers match)                    ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  VARIANTS FROM rCRS                                                     |
|  +--------------+------------+------------+-------------------------+   |
|  | Position     | rCRS       | Sample     | Region                  |   |
|  +--------------+------------+------------+-------------------------+   |
|  | 263          | A          | G          | HVS2                    |   |
|  | 315.1        | -          | C          | HVS2                    |   |
|  | 750          | A          | G          | Coding                  |   |
|  | ...          |            |            |                         |   |
|  +--------------+------------+------------+-------------------------+   |
|                                                                         |
|  [Export FASTA]  [View Full Sequence]                                   |
|                                                                         |
+-------------------------------------------------------------------------+
```

### Subject Tab: Ancestry Composition

```
+-------------------------------------------------------------------------+
|  ANCESTRY COMPOSITION                                   [Run Analysis]  |
+-------------------------------------------------------------------------+
|                                                                         |
|  +---------------------------------------------------------------------+|
|  |                                                                     ||
|  |   Northwestern European  ====================     78.2%             ||
|  |   +- British & Irish     ================         54.1%             ||
|  |   +- French & German     ======                   16.8%             ||
|  |   +- Scandinavian        ===                       7.3%             ||
|  |                                                                     ||
|  |   Southern European      =====                    12.4%             ||
|  |   +- Italian             ====                      9.1%             ||
|  |   +- Iberian             =                         3.3%             ||
|  |                                                                     ||
|  |   Eastern European       ===                       6.2%             ||
|  |                                                                     ||
|  |   Broadly European       =                         3.2%             ||
|  |                                                                     ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  Reference Panel: v5.2 (Dec 2024)  |  Confidence: 90%                   |
|                                                                         |
+-------------------------------------------------------------------------+
```

### Subject Tab: IBD Matches

All IBD matches for this subject:

```
+-------------------------------------------------------------------------+
|  IBD MATCHES                                 [Run Match] [Import Matches]|
+-------------------------------------------------------------------------+
|                                                                         |
|  [Filter matches...]             [Min cM: 20 v]    [Show: All v]        |
|                                                                         |
|  +---------------------------------------------------------------------+|
|  | Match         | Shared   | Segments | Longest | Relationship        ||
|  +---------------+----------+----------+---------+---------------------+|
|  | * Bob Smith   | 1,847 cM | 42       | 287 cM  | Close Family        ||
|  |   Jane Doe    | 847 cM   | 28       | 94 cM   | 1st-2nd Cousin      ||
|  |   Tom Wilson  | 127 cM   | 8        | 32 cM   | 3rd-4th Cousin      ||
|  |   Sue Brown   | 45 cM    | 3        | 18 cM   | 4th-5th Cousin      ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  SELECTED: Bob Smith                                                    |
|  +---------------------------------------------------------------------+|
|  | CHROMOSOME BROWSER                                                  ||
|  | Chr 1  [====    ====                              ========       ]  ||
|  | Chr 2  [    ========              ====                           ]  ||
|  | Chr 3  [                ============                             ]  ||
|  | ...                                                                 ||
|  | [=] John Smith (this subject)    [#] Bob Smith                      ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  [View in Compare Mode]  [Export Segment Data]                          |
|                                                                         |
+-------------------------------------------------------------------------+
```

### Subject Tab: Data Sources

Raw data management:

```
+-------------------------------------------------------------------------+
|  DATA SOURCES                                               [+ Add Data]|
+-------------------------------------------------------------------------+
|                                                                         |
|  SEQUENCING RUNS                                                        |
|  +---------------------------------------------------------------------+|
|  | v WGS_NovaSeq_2024 (Primary)                               [...]    ||
|  |   Platform: Illumina NovaSeq 6000                                   ||
|  |   +- hg38 alignment - 32.4x mean coverage - GATK 4.6               ||
|  |      Analysis: [x] Metrics  [x] Callable  [x] Y-DNA  [x] mtDNA     ||
|  |                                                                     ||
|  | > WGS_HiSeq_2022                                           [...]    ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  CHIP / ARRAY PROFILES                                                  |
|  +---------------------------------------------------------------------+|
|  |   23andMe v5 Import (Dec 2024)                             [...]    ||
|  |   SNPs: 642,824 - Call rate: 99.2%                                  ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  STR PROFILES                                                           |
|  +---------------------------------------------------------------------+|
|  |   Y-111 Panel (FTDNA Import)                               [...]    ||
|  |   111 markers - DYS393: 13, DYS390: 24, ...                         ||
|  +---------------------------------------------------------------------+|
|                                                                         |
+-------------------------------------------------------------------------+
```

---

## Tab 3: Projects

Project management with member lists and drag-drop support.

```
+-------------------------------------------------------------------------+
| PROJECTS                                                  [+ New Project]|
+-------------------------------------------------------------------------+
| [Search projects...]                                                    |
+-------------------------------------------------------------------------+
|                                                                         |
|  +---------------------------------------------------------------------+|
|  | Smith Family Research                                    [Edit] [x] ||
|  | Created: Nov 2024 | Members: 12 | Admin: jkane                      ||
|  +---------------------------------------------------------------------+|
|  | Doe Genetic Study                                        [Edit] [x] ||
|  | Created: Oct 2024 | Members: 8 | Admin: jkane                       ||
|  +---------------------------------------------------------------------+|
|  | Asia Pacific Project                                     [Edit] [x] ||
|  | Created: Dec 2024 | Members: 5 | Admin: jkane                       ||
|  +---------------------------------------------------------------------+|
|                                                                         |
+-------------------------------------------------------------------------+

SELECTED PROJECT: Smith Family Research

+-------------------------------------------------------------------------+
|  PROJECT MEMBERS                                         [+ Add Member] |
+-------------------------------------------------------------------------+
|  +---------------------------------------------------------------------+|
|  | John Smith    | DU-001 | R1b-P312 | H1a1  | [View] [Remove]         ||
|  | Mary Smith    | DU-006 | -        | H1a1  | [View] [Remove]         ||
|  | Bob Smith     | DU-002 | R1b-P312 | H1a1  | [View] [Remove]         ||
|  | ...                                                                 ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  +- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+|
|  |              Drag subjects here to add to project                   ||
|  +- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+|
|                                                                         |
+-------------------------------------------------------------------------+
```

---

## Compare Mode

Accessed via multi-select in Subject grid or "Compare with Others" button:

```
+-------------------------------------------------------------------------+
|  COMPARE SUBJECTS                                          [+ Add] [x]  |
+---------------------+---------------------+-----------------------------+
|  John Smith (DU-001)|  Bob Smith (DU-002) |  [Empty - drag to add]      |
+---------------------+---------------------+-----------------------------+
|  Comparison Type: [Y-DNA v]  [STR v]  [Ancestry v]  [IBD v]             |
+-------------------------------------------------------------------------+
|                                                                         |
|  Y-STR COMPARISON (67 markers)                                          |
|  +---------------------------------------------------------------------+|
|  | Marker   | John Smith | Bob Smith  | Status                         ||
|  +----------+------------+------------+--------------------------------+|
|  | DYS393   | 13         | 13         | [x] Match                      ||
|  | DYS390   | 24         | 24         | [x] Match                      ||
|  | DYS19    | 14         | 15         | [!] Diff (1 step)              ||
|  | DYS391   | 11         | 11         | [x] Match                      ||
|  | ...      |            |            |                                ||
|  +---------------------------------------------------------------------+|
|                                                                         |
|  SUMMARY: Genetic Distance = 1 (out of 67 markers)                      |
|  Estimated MRCA: 4-6 generations                                        |
|                                                                         |
+-------------------------------------------------------------------------+
```

---

## Design System

### Color Palette (Modern Dark Theme)

```
Background Hierarchy:
  --bg-primary:    #1a1a2e    (main canvas)
  --bg-secondary:  #16213e    (cards/panels)
  --bg-tertiary:   #0f3460    (elevated elements)

Accent Colors:
  --accent-y-dna:  #4ade80    (green - Y chromosome)
  --accent-mt-dna: #60a5fa    (blue - mitochondrial)
  --accent-ibd:    #f472b6    (pink - IBD matching)
  --accent-warn:   #fbbf24    (amber - warnings)
  --accent-error:  #f87171    (red - errors)

Text:
  --text-primary:  #f1f5f9
  --text-secondary:#94a3b8
  --text-muted:    #64748b
```

### Typography

```
Headers:     Inter or SF Pro Display (system fonts)
Body:        Inter or system-ui
Monospace:   JetBrains Mono (for IDs, sequences)
```

### Component Guidelines

1. **Cards with subtle gradients** instead of flat boxes
2. **Rounded corners** (8-12px radius)
3. **Subtle shadows** for depth hierarchy
4. **Icon integration** - genetic/science iconography
5. **Progress indicators** - ring spinners, skeleton loaders
6. **Toast notifications** - non-blocking feedback

---

## Implementation Phases

### Phase 1: Foundation
1. Implement tab-based navigation shell
2. Create new dashboard home view
3. Upgrade data grid with multi-select and column customization
4. Add batch operation support

### Phase 2: Enhanced Subject Detail
1. Refactor subject detail into tabbed interface
2. Move analysis results into dedicated sub-tabs
3. Add IBD matches tab
4. Implement ancestry composition view

### Phase 3: Comparison & Visualization
1. Haplogroup tree visualization (interactive SVG/Canvas)
2. Comparison view for multiple subjects
3. Chromosome browser component
4. Enhanced STR comparison tools

### Phase 4: IBD & Ancestry (when features ready)
1. IBD matching integration
2. Relationship estimation display
3. Ancestry composition visualization
4. Cross-subject IBD network view

---

## Key Architectural Changes

| Current | Proposed |
|---------|----------|
| Single SplitPane | TabPane with specialized views |
| ListView for subjects | TableView with virtual scrolling |
| Modal dialogs for results | Inline expandable panels + sub-tabs |
| Text-only haplogroups | Interactive tree visualization |
| One detail view | Comparison view + rich detail |
| Imperative styling | CSS variables + theme system |
| 5 top-level concepts | 3 entity-focused tabs |
| Analysis as global feature | Analysis within Subject detail |
| Matching as global feature | IBD Matches within Subject detail |
