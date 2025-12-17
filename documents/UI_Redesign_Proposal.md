# UI Redesign Proposal: DUNavigator

**Date:** December 2024
**Last Updated:** December 17, 2024
**Status:** In Progress - Phases 1-3 Complete
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

## Implementation Status

### Current Implementation (December 2024)

#### Completed Features

**Navigation & Layout**
- [x] Tab-based navigation shell (Dashboard, Subjects, Projects)
- [x] Modern dark theme with color palette
- [x] Status bar with workspace info
- [x] Responsive layout with proper spacing

**Dashboard Tab**
- [x] Workspace summary cards (subject count, project count, Y-DNA count, mtDNA count, IBD count)
- [x] Pending work queue with clickable items
- [x] Haplogroup distribution charts (Y-DNA and mtDNA)
- [x] Recent activity feed

**Subjects Tab**
- [x] Data grid with multi-select support
- [x] Column customization (show/hide columns)
- [x] Search and filtering
- [x] Batch operations toolbar (Compare, Batch Analyze, Add to Project)
- [x] Inline status indicators with color coding

**Subject Detail View**
- [x] Tabbed interface (Overview, Y-DNA, mtDNA, Ancestry, IBD, Data Sources)
- [x] Header with subject info and action buttons

**Overview Tab**
- [x] Y-DNA haplogroup card with confidence and quality badge
- [x] mtDNA haplogroup card with confidence and quality badge
- [x] Y-STR summary panel with FTDNA panel indicators (Y-12 through Y-700)
- [x] Ancestry placeholder card
- [x] IBD matches placeholder card
- [x] Data sources summary (sequencing runs, chip profiles, STR profiles)

**Y-DNA Tab**
- [x] Terminal haplogroup display with phylogenetic path
- [x] Derived/Ancestral SNP counts
- [x] Confidence level with color coding
- [x] Quality rating with color coding
- [x] Last analyzed timestamp
- [x] Source reconciliation table (multi-run comparison)
- [x] Y-Chromosome ideogram with region bands (PAR, X-degenerate, Ampliconic, etc.)
- [x] View Full Y Profile button (opens detailed dialog with variant markers)
- [x] Run Analysis button with slide-in progress panel

**mtDNA Tab**
- [x] Haplogroup display with phylogenetic path
- [x] Confidence level with color coding
- [x] Quality rating with color coding
- [x] Last analyzed timestamp
- [x] Source reconciliation table (multi-run comparison)
- [x] Variants from rCRS table with region classification (HVS1/HVS2/HVS3/Coding)
- [x] Export FASTA placeholder button
- [x] Run Analysis button

**Data Sources Tab**
- [x] Sequencing runs list with expandable details
- [x] Chip/Array profiles list
- [x] STR profiles list
- [x] Add Data button with file type auto-detection
- [x] VCF import with metadata dialog

**Compare View**
- [x] Subject comparison dialog (2-3 subjects)
- [x] Y-STR comparison tab with genetic distance calculation
- [x] Y-DNA haplogroup comparison tab
- [x] mtDNA haplogroup comparison tab
- [x] Match/mismatch highlighting with step difference display

**Projects Tab**
- [x] Project list with member counts
- [x] Project detail with member grid
- [x] Add/remove members
- [x] Create new project

**Analysis Features**
- [x] Slide-in progress panel for running analyses
- [x] Progress tracking with task name and status
- [x] Cancel button for running analyses

**Internationalization**
- [x] Full i18n support with messages.properties
- [x] Locale-aware number/date formatting via Formatters

#### Pending Features (Future Work)

**Ancestry Tab** (Phase 4)
- [ ] Ancestry composition visualization (bar chart)
- [ ] Reference panel information
- [ ] Run Analysis integration (requires backend implementation)

**IBD Matches Tab** (Phase 4)
- [ ] IBD matches table with filtering (UI placeholder exists)
- [ ] Chromosome browser visualization
- [ ] Relationship estimation display
- [ ] Import Matches functionality (requires backend implementation)
- [ ] Run Match functionality (requires backend implementation)

**Advanced Visualizations** (Future)
- [ ] Interactive haplogroup tree visualization (SVG/Canvas)
- [x] Y-Chromosome ideogram with region bands (moved from Y Profile dialog to Y-DNA tab)
- [ ] Y-Chromosome ideogram with SNP position markers (available in View Full Y Profile)
- [ ] Cross-subject IBD network view

---

## Implementation Phases

### Phase 1: Foundation - COMPLETE
1. [x] Implement tab-based navigation shell
2. [x] Create new dashboard home view
3. [x] Upgrade data grid with multi-select and column customization
4. [x] Add batch operation support

### Phase 2: Enhanced Subject Detail - COMPLETE
1. [x] Refactor subject detail into tabbed interface
2. [x] Move analysis results into dedicated sub-tabs
3. [x] Add IBD matches tab (placeholder UI)
4. [x] Implement ancestry composition view (placeholder UI)

### Phase 3: Comparison & Visualization - COMPLETE
1. [ ] Haplogroup tree visualization (interactive SVG/Canvas) - DEFERRED
2. [x] Comparison view for multiple subjects
3. [ ] Chromosome browser component - DEFERRED (awaiting IBD backend)
4. [x] Enhanced STR comparison tools

### Phase 4: IBD & Ancestry - FUTURE WORK
*Requires backend implementation before UI can be fully functional*

1. [ ] IBD matching integration
   - Backend: IBD segment detection algorithm
   - Backend: Match storage and retrieval
   - UI: Wire up matches table
   - UI: Implement chromosome browser

2. [ ] Relationship estimation display
   - Backend: cM to relationship mapping
   - UI: Relationship badges in matches table

3. [ ] Ancestry composition visualization
   - Backend: Ancestry inference algorithm
   - Backend: Reference panel integration
   - UI: Stacked bar chart component
   - UI: Confidence intervals display

4. [ ] Cross-subject IBD network view
   - Backend: Network graph data structure
   - UI: Force-directed graph visualization

---

## Key Architectural Changes

| Original Design | V2 Implementation | Status |
|-----------------|-------------------|--------|
| Single SplitPane | TabPane with specialized views | DONE |
| ListView for subjects | TableView with virtual scrolling | DONE |
| Modal dialogs for results | Inline expandable panels + sub-tabs | DONE |
| Text-only haplogroups | Interactive tree visualization | DEFERRED |
| One detail view | Comparison view + rich detail | DONE |
| Imperative styling | CSS variables + theme system | DONE |
| 5 top-level concepts | 3 entity-focused tabs | DONE |
| Analysis as global feature | Analysis within Subject detail | DONE |
| Matching as global feature | IBD Matches within Subject detail | PLACEHOLDER |

---

## Planned UI Enhancements

### Near-Term Improvements

These improvements can be added without backend changes:

**Y-DNA Tab Enhancements**
- [ ] SNP list export button (CSV/TSV format)
- [ ] Haplogroup badge with phylogenetic depth indicator
- [ ] Tree provider selector (switch between FTDNA/DecodingUs trees)
- [ ] Private SNP discovery panel (when available from analysis)

**mtDNA Tab Enhancements**
- [ ] Heteroplasmy indicator for variants
- [ ] Region highlighting in variant table (color-coded by region)
- [ ] Export FASTA implementation (construct sequence from variants)
- [ ] mtDNA sequence viewer (circular or linear visualization)

**Compare View Enhancements**
- [ ] Export comparison report (PDF/HTML)
- [ ] MRCA estimation based on genetic distance
- [ ] Modal STR difference highlighting
- [ ] Add more than 3 subjects (scrollable comparison)

**Dashboard Enhancements**
- [ ] Clickable haplogroup distribution bars (filter subjects)
- [ ] Analysis queue progress indicators
- [ ] Quick actions (add subject, import data)

**Data Sources Enhancements**
- [ ] Data quality summary per source
- [ ] Coverage histogram visualization
- [ ] Batch import from folder

### Future UI Work (Requires Backend)

**Ancestry Composition** (Phase 4)
```
+---------------------------------------------------------------------+
|  ANCESTRY COMPOSITION                                               |
+---------------------------------------------------------------------+
|                                                                     |
|  [Stacked horizontal bar chart]                                     |
|  ================================================================== |
|  | NW European 78% | S. European 12% | E. European 6% | Other 4% | |
|  ================================================================== |
|                                                                     |
|  DETAILED BREAKDOWN                              Confidence: 90%    |
|  +----------------------------+                                     |
|  | British & Irish     54.1% | ===============                      |
|  | French & German     16.8% | =====                                |
|  | Scandinavian         7.3% | ==                                   |
|  | Italian              9.1% | ===                                  |
|  | Iberian              3.3% | =                                    |
|  | Eastern European     6.2% | ==                                   |
|  | Broadly European     3.2% | =                                    |
|  +----------------------------+                                     |
|                                                                     |
|  Reference Panel: v5.2 (Dec 2024)                                   |
|  [Run Analysis]  [Export Report]                                    |
+---------------------------------------------------------------------+
```

**IBD Matches** (Phase 4)
```
+---------------------------------------------------------------------+
|  IBD MATCHES                            [Run Match] [Import Matches]|
+---------------------------------------------------------------------+
|                                                                     |
|  FILTERS                                                            |
|  Min cM: [====o====] 20 cM     Relationship: [All v]                |
|  Show: [x] In workspace  [ ] External matches                       |
|                                                                     |
|  MATCHES (47 found)                                                 |
|  +---------------------+--------+------+-------+------------------+ |
|  | Match               | Shared | Segs | Long. | Relationship     | |
|  +---------------------+--------+------+-------+------------------+ |
|  | * Bob Smith [ws]    | 1,847  | 42   | 287   | Close Family     | |
|  |   Jane Doe [ws]     | 847    | 28   | 94    | 1st-2nd Cousin   | |
|  |   External Match 1  | 127    | 8    | 32    | 3rd-4th Cousin   | |
|  +---------------------+--------+------+-------+------------------+ |
|  [ws] = In workspace                                                |
|                                                                     |
|  CHROMOSOME BROWSER                                                 |
|  +-------------------------------------------------------------------+
|  | Chr 1  [====    ====                              ========     ] |
|  | Chr 2  [    ========              ====                         ] |
|  | Chr 3  [                ============                           ] |
|  | Chr 4  [  ====                                    ====         ] |
|  | ...                                                              |
|  +-------------------------------------------------------------------+
|  Legend: [===] Shared IBD segment                                   |
|                                                                     |
|  [View in Compare Mode]  [Export Segments]  [Add to Project]        |
+---------------------------------------------------------------------+
```

**Haplogroup Tree Visualization** (Future)
```
+---------------------------------------------------------------------+
|  Y-DNA PHYLOGENETIC TREE                           [Zoom] [Reset]   |
+---------------------------------------------------------------------+
|                                                                     |
|  R ─┬─ R1 ─┬─ R1a ─── ...                                          |
|     │      │                                                        |
|     │      └─ R1b ─┬─ R1b-M269 ─┬─ R1b-L151 ─┬─ R1b-P312 ●         |
|     │              │            │            │                      |
|     │              │            │            └─ R1b-U106 ───...     |
|     │              │            │                                   |
|     │              │            └─ R1b-L11 ───...                   |
|     │              │                                                |
|     │              └─ R1b-V88 ───...                                |
|     │                                                               |
|     └─ R2 ───...                                                    |
|                                                                     |
|  ● = Current subject position                                       |
|  Click any node to see details                                      |
|                                                                     |
+---------------------------------------------------------------------+
```
