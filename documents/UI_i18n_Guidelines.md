# Internationalization (i18n) Guidelines

**Date:** December 2024
**Status:** Draft Specification
**Related:** UI_Redesign_Proposal.md

---

## Overview

For a genetic genealogy application with global reach, i18n must be a first-class concern. This document outlines the architecture, implementation patterns, and guidelines for supporting multiple languages in DUNavigator.

---

## Target Locales

### Priority 1 (Initial Release)
- **English** (en) - Default
- **German** (de) - Strong genealogy community
- **Spanish** (es) - Latin America, Spain

### Priority 2
- **French** (fr) - Canada, France
- **Portuguese** (pt-BR) - Brazil

### Priority 3
- **Japanese** (ja)
- **Simplified Chinese** (zh-CN)

### Priority 4 (RTL Support)
- **Arabic** (ar)
- **Hebrew** (he)

---

## Resource Bundle Structure

```
src/main/resources/
+-- i18n/
|   +-- messages.properties          # Default (English)
|   +-- messages_de.properties       # German
|   +-- messages_es.properties       # Spanish
|   +-- messages_fr.properties       # French
|   +-- messages_pt_BR.properties    # Brazilian Portuguese
|   +-- messages_ja.properties       # Japanese
|   +-- messages_zh_CN.properties    # Simplified Chinese
|   +-- messages_ar.properties       # Arabic (RTL)
|   +-- messages_he.properties       # Hebrew (RTL)
|   +-- haplogroups.properties       # Haplogroup descriptions (shared)
|   +-- scientific_terms.properties  # Domain terminology (shared)
+-- style/
    +-- style.css                    # Base styles
    +-- style-rtl.css                # RTL overrides
```

---

## Message Key Convention

Use hierarchical dot-notation keys organized by feature area.

### Key Naming Rules

1. Use lowercase with dots as separators
2. Group by feature/component first
3. Use descriptive, semantic names
4. Keep keys reasonably short

### Message File Template

```properties
# =============================================================================
# messages.properties - English (Default)
# =============================================================================

# -----------------------------------------------------------------------------
# Navigation
# -----------------------------------------------------------------------------
nav.dashboard=Dashboard
nav.subjects=Subjects
nav.projects=Projects

# -----------------------------------------------------------------------------
# Dashboard
# -----------------------------------------------------------------------------
dashboard.title=Workspace Overview
dashboard.subjects.count={0} Subjects
dashboard.subjects.with_ydna={0} with Y-DNA
dashboard.subjects.with_mtdna={0} with mtDNA
dashboard.subjects.with_ibd={0} with IBD data
dashboard.projects.count={0} Projects
dashboard.pending_work=Pending Work
dashboard.pending_analyses={0} Pending Analyses
dashboard.run_all=Run All
dashboard.recent_activity=Recent Activity
dashboard.haplogroup_distribution=Haplogroup Distributions
dashboard.ydna_distribution=Y-DNA Distribution
dashboard.mtdna_distribution=mtDNA Distribution

# -----------------------------------------------------------------------------
# Subject Grid
# -----------------------------------------------------------------------------
subjects.title=Subjects
subjects.add=Add New Subject
subjects.search.placeholder=Search name, ID, haplogroup...
subjects.filter=Filters
subjects.columns=Columns
subjects.selected={0} selected
subjects.compare=Compare
subjects.batch_analyze=Batch Analyze
subjects.add_to_project=Add to Project

# -----------------------------------------------------------------------------
# Column Headers
# -----------------------------------------------------------------------------
column.id=ID
column.name=Name
column.ydna=Y-DNA
column.mtdna=mtDNA
column.project=Project
column.status=Status
column.sex=Sex
column.date_added=Added

# -----------------------------------------------------------------------------
# Subject Detail Tabs
# -----------------------------------------------------------------------------
subject.tab.overview=Overview
subject.tab.ydna=Y-DNA
subject.tab.mtdna=mtDNA
subject.tab.ancestry=Ancestry
subject.tab.ibd=IBD Matches
subject.tab.data=Data Sources

# -----------------------------------------------------------------------------
# Status Indicators
# -----------------------------------------------------------------------------
status.complete=Complete
status.pending=Pending
status.error=Error
status.none=None
status.consistent=Consistent
status.inconsistent=Inconsistent

# -----------------------------------------------------------------------------
# Haplogroup Panel
# -----------------------------------------------------------------------------
haplogroup.title=Haplogroup
haplogroup.ydna.title=Y-DNA
haplogroup.mtdna.title=mtDNA
haplogroup.derived=Derived
haplogroup.ancestral=Ancestral
haplogroup.callable=Callable
haplogroup.confidence=Confidence
haplogroup.confidence.high=HIGH
haplogroup.confidence.medium=MEDIUM
haplogroup.confidence.low=LOW
haplogroup.view_profile=View Full Y Profile
haplogroup.view_details=View mtDNA Details
haplogroup.terminal=Terminal Haplogroup
haplogroup.phylogenetic_path=Phylogenetic Path
haplogroup.confirmed=Haplogroup confirmed

# -----------------------------------------------------------------------------
# Analysis
# -----------------------------------------------------------------------------
analysis.title=Analysis
analysis.not_analyzed=Not yet analyzed
analysis.run=Run Analysis
analysis.pending=Analysis pending
analysis.complete=Analysis complete
analysis.failed=Analysis failed
analysis.retry=Retry Analysis
analysis.coverage=Coverage
analysis.callable_loci=Callable Loci
analysis.wgs_metrics=WGS Metrics
analysis.source_reconciliation=Source Reconciliation
analysis.quality=Quality
analysis.quality.excellent=Excellent
analysis.quality.good=Good
analysis.quality.fair=Fair
analysis.quality.poor=Poor

# -----------------------------------------------------------------------------
# Ancestry Composition
# -----------------------------------------------------------------------------
ancestry.title=Ancestry Composition
ancestry.reference_panel=Reference Panel
ancestry.confidence=Confidence

# -----------------------------------------------------------------------------
# IBD Matching
# -----------------------------------------------------------------------------
ibd.title=IBD Matching
ibd.matches=IBD Matches
ibd.run_match=Run Match
ibd.import_matches=Import Matches
ibd.filter_matches=Filter matches...
ibd.min_cm=Min cM
ibd.show_all=All
ibd.matches_found={0} matches found
ibd.no_matches=No matches found
ibd.closest_match=Closest: {0} ({1} cM)
ibd.shared_cm=Shared cM
ibd.segments=Segments
ibd.longest=Longest
ibd.relationship=Relationship Est.
ibd.chromosome_browser=Chromosome Browser
ibd.export_segments=Export Segment Data
ibd.view_compare=View in Compare Mode

# -----------------------------------------------------------------------------
# Relationship Estimates
# -----------------------------------------------------------------------------
relationship.close_family=Close Family
relationship.first_cousin=1st Cousin
relationship.second_cousin=2nd Cousin
relationship.first_second_cousin=1st-2nd Cousin
relationship.third_cousin=3rd Cousin
relationship.third_fourth_cousin=3rd-4th Cousin
relationship.fourth_cousin=4th Cousin
relationship.fourth_fifth_cousin=4th-5th Cousin
relationship.distant=Distant Relative

# -----------------------------------------------------------------------------
# Data Sources
# -----------------------------------------------------------------------------
data.title=Data Sources
data.add=Add Data
data.sequencing_runs=Sequencing Runs
data.chip_profiles=Chip / Array Profiles
data.str_profiles=STR Profiles
data.primary=Primary
data.platform=Platform
data.alignment=Alignment
data.coverage=Coverage
data.snps=SNPs
data.call_rate=Call rate
data.markers=Markers

# -----------------------------------------------------------------------------
# Compare Mode
# -----------------------------------------------------------------------------
compare.title=Compare Subjects
compare.add=Add Subject
compare.empty=Empty - drag to add
compare.type.ydna=Y-DNA
compare.type.str=STR
compare.type.ancestry=Ancestry
compare.type.ibd=IBD
compare.genetic_distance=Genetic Distance
compare.mrca_estimate=Estimated MRCA
compare.match=Match
compare.diff=Diff
compare.step=step
compare.steps=steps
compare.summary=Summary
compare.out_of=out of
compare.generations=generations

# -----------------------------------------------------------------------------
# Projects
# -----------------------------------------------------------------------------
projects.title=Projects
projects.add=New Project
projects.search=Search projects...
projects.members=Members
projects.admin=Admin
projects.created=Created
projects.add_member=Add Member
projects.remove_member=Remove
projects.view_member=View
projects.drag_hint=Drag subjects here to add to project

# -----------------------------------------------------------------------------
# Actions
# -----------------------------------------------------------------------------
action.edit=Edit
action.delete=Delete
action.save=Save
action.cancel=Cancel
action.close=Close
action.back=Back to List
action.analyze=Analyze
action.export=Export
action.import=Import
action.view=View
action.remove=Remove
action.refresh=Refresh

# -----------------------------------------------------------------------------
# Confirmations
# -----------------------------------------------------------------------------
confirm.delete.title=Confirm Delete
confirm.delete.subject=Are you sure you want to delete {0}?
confirm.delete.project=Delete project "{0}" and remove all member associations?
confirm.unsaved=You have unsaved changes. Discard them?

# -----------------------------------------------------------------------------
# Errors
# -----------------------------------------------------------------------------
error.generic=An error occurred
error.load_failed=Failed to load {0}
error.save_failed=Failed to save {0}
error.analysis_failed=Analysis failed: {0}
error.connection=Connection error. Working offline.
error.not_found={0} not found
error.invalid_input=Invalid input: {0}

# -----------------------------------------------------------------------------
# Time and Date
# -----------------------------------------------------------------------------
date.added=Added {0}
date.updated=Updated {0}
date.created=Created {0}
time.ago.hours={0}h ago
time.ago.days={0}d ago
time.yesterday=Yesterday
time.today=Today

# -----------------------------------------------------------------------------
# Units (symbols are universal, labels may be translated)
# -----------------------------------------------------------------------------
unit.centimorgans=cM
unit.coverage=x
unit.megabases=Mb
unit.percent=%
```

---

## ScalaFX i18n Implementation

### Core I18n Object

```scala
// src/main/scala/com/decodingus/i18n/I18n.scala

package com.decodingus.i18n

import java.util.{Locale, ResourceBundle, MissingResourceException}
import java.text.MessageFormat
import scalafx.beans.property.{ObjectProperty, StringProperty, ReadOnlyStringProperty}
import scalafx.beans.binding.Bindings

object I18n {

  // Observable locale for reactive UI updates
  val currentLocale: ObjectProperty[Locale] = ObjectProperty(Locale.getDefault)

  // Cached bundle reference
  private var cachedBundle: ResourceBundle = loadBundle()

  private def loadBundle(): ResourceBundle = {
    ResourceBundle.getBundle("i18n.messages", currentLocale.value)
  }

  // Reload bundle when locale changes
  currentLocale.onChange { (_, _, _) =>
    cachedBundle = loadBundle()
  }

  /**
   * Get a simple translated string.
   * Returns "!key!" if key is missing (makes missing translations visible).
   */
  def t(key: String): String = {
    try {
      cachedBundle.getString(key)
    } catch {
      case _: MissingResourceException =>
        System.err.println(s"Missing i18n key: $key")
        s"!$key!"
    }
  }

  /**
   * Get a translated string with parameters (ICU MessageFormat).
   * Example: t("dashboard.subjects.count", 47) -> "47 Subjects"
   */
  def t(key: String, args: Any*): String = {
    val pattern = t(key)
    if (args.isEmpty) pattern
    else MessageFormat.format(pattern, args.map(_.asInstanceOf[AnyRef]): _*)
  }

  /**
   * Create a reactive StringProperty that updates when locale changes.
   * Use for static labels that don't have parameters.
   */
  def bind(key: String): StringProperty = {
    val prop = StringProperty(t(key))
    currentLocale.onChange { (_, _, _) =>
      prop.value = t(key)
    }
    prop
  }

  /**
   * Create a reactive binding with dynamic parameters.
   * The args function is re-evaluated on each locale change.
   */
  def bind(key: String, args: () => Seq[Any]): StringProperty = {
    val prop = StringProperty(t(key, args(): _*))
    currentLocale.onChange { (_, _, _) =>
      prop.value = t(key, args(): _*)
    }
    prop
  }

  /**
   * Create a read-only binding for use in UI components.
   */
  def bindReadOnly(key: String): ReadOnlyStringProperty = bind(key)

  /**
   * Check if current locale is RTL (right-to-left).
   */
  def isRTL: Boolean = {
    val rtlLanguages = Set("ar", "he", "fa", "ur")
    rtlLanguages.contains(currentLocale.value.getLanguage)
  }

  /**
   * Switch locale at runtime.
   */
  def setLocale(locale: Locale): Unit = {
    currentLocale.value = locale
    Locale.setDefault(locale)
  }

  /**
   * Get list of supported locales for UI picker.
   */
  def supportedLocales: Seq[Locale] = Seq(
    Locale.ENGLISH,
    Locale.GERMAN,
    new Locale("es"),
    Locale.FRENCH,
    new Locale("pt", "BR"),
    Locale.JAPANESE,
    Locale.SIMPLIFIED_CHINESE,
    new Locale("ar"),
    new Locale("he")
  )
}
```

### Locale-Aware Formatters

```scala
// src/main/scala/com/decodingus/i18n/Formatters.scala

package com.decodingus.i18n

import java.text.NumberFormat
import java.time.LocalDate
import java.time.format.{DateTimeFormatter, FormatStyle}

object Formatters {

  /**
   * Format a number with locale-appropriate separators.
   * Example: 1234567 -> "1,234,567" (en) or "1.234.567" (de)
   */
  def formatNumber(n: Number): String = {
    NumberFormat.getNumberInstance(I18n.currentLocale.value).format(n)
  }

  /**
   * Format a decimal number with specified precision.
   */
  def formatDecimal(d: Double, decimals: Int = 1): String = {
    val formatter = NumberFormat.getNumberInstance(I18n.currentLocale.value)
    formatter.setMinimumFractionDigits(decimals)
    formatter.setMaximumFractionDigits(decimals)
    formatter.format(d)
  }

  /**
   * Format a percentage.
   * Example: 0.982 -> "98.2%" (en) or "98,2 %" (fr)
   */
  def formatPercent(d: Double): String = {
    NumberFormat.getPercentInstance(I18n.currentLocale.value).format(d)
  }

  /**
   * Format a date in medium style.
   * Example: "Dec 15, 2024" (en) or "15. Dez. 2024" (de)
   */
  def formatDate(date: LocalDate): String = {
    DateTimeFormatter
      .ofLocalizedDate(FormatStyle.MEDIUM)
      .withLocale(I18n.currentLocale.value)
      .format(date)
  }

  /**
   * Format a date in short style.
   */
  def formatDateShort(date: LocalDate): String = {
    DateTimeFormatter
      .ofLocalizedDate(FormatStyle.SHORT)
      .withLocale(I18n.currentLocale.value)
      .format(date)
  }

  /**
   * Format coverage value with 'x' suffix.
   * Coverage symbol is universal, but number formatting is locale-specific.
   */
  def formatCoverage(coverage: Double): String = {
    s"${formatDecimal(coverage)}x"
  }

  /**
   * Format centimorgan value.
   * cM symbol is universal.
   */
  def formatCentimorgans(cm: Double): String = {
    s"${formatNumber(cm.toLong)} cM"
  }
}
```

---

## Usage Patterns

### Static Labels

```scala
import com.decodingus.i18n.I18n.{t, bind}

// One-time evaluation (doesn't update on locale change)
val deleteButton = new Button(t("action.delete"))

// Reactive binding (updates when locale changes)
val titleLabel = new Label {
  text <== bind("dashboard.title")
}
```

### Parameterized Messages

```scala
// One-time with parameters
val message = t("dashboard.subjects.count", subjects.size)
// -> "47 Subjects"

// Reactive with dynamic parameters
val countLabel = new Label {
  text <== bind("dashboard.subjects.count", () => Seq(viewModel.subjects.size))
}

// Update manually when data changes
viewModel.subjects.onChange { (_, _) =>
  countLabel.text = t("dashboard.subjects.count", viewModel.subjects.size)
}
```

### Tab Labels with Icons

```scala
private def createTab(key: String, icon: String, content: Node): Tab = new Tab {
  // Combine icon with translated text
  text <== bind(key).map(s => s"$icon $s")
  closable = false
  this.content = content
}

val tabPane = new TabPane {
  tabs = Seq(
    createTab("nav.dashboard", "\uD83D\uDCCA", dashboardContent),  // chart icon
    createTab("nav.subjects", "\uD83D\uDC65", subjectsContent),    // people icon
    createTab("nav.projects", "\uD83D\uDCC1", projectsContent)     // folder icon
  )
}
```

### Table Columns

```scala
val nameColumn = new TableColumn[Subject, String] {
  text <== bind("column.name")
  cellValueFactory = { p => StringProperty(p.value.name) }
}

val statusColumn = new TableColumn[Subject, String] {
  text <== bind("column.status")
  cellValueFactory = { p =>
    // Translate status values
    StringProperty(t(s"status.${p.value.status.toLowerCase}"))
  }
}
```

---

## UI Layout Guidelines for i18n

### Text Expansion

German and French text is typically **30-40% longer** than English:

| Language | "Subjects" | "Analysis" | "Run All" |
|----------|------------|------------|-----------|
| English  | Subjects   | Analysis   | Run All   |
| German   | Probanden  | Analyse    | Alle ausführen |
| French   | Sujets     | Analyse    | Tout exécuter |
| Spanish  | Sujetos    | Análisis   | Ejecutar todo |

### Layout Rules

1. **Never use fixed-width containers for text**

```scala
// BAD - Fixed width will clip German text
val label = new Label(t("nav.subjects")) {
  prefWidth = 80
}

// GOOD - Flexible width with minimum
val label = new Label {
  text <== bind("nav.subjects")
  minWidth = Region.USE_PREF_SIZE
  maxWidth = Double.MaxValue
}
```

2. **Use HBox/VBox with growth priorities**

```scala
val buttonBar = new HBox(10) {
  children = Seq(
    new Button { text <== bind("action.save") },
    new Region { hgrow = Priority.Always }, // Spacer
    new Button { text <== bind("action.cancel") }
  )
}
```

3. **Labels above fields (not beside) for forms**

```scala
// Better for text expansion
val formField = new VBox(5) {
  children = Seq(
    new Label { text <== bind("field.email") },
    new TextField()
  )
}
```

4. **Allow buttons to grow**

```scala
val button = new Button {
  text <== bind("action.batch_analyze")
  minWidth = Region.USE_PREF_SIZE
  maxWidth = Region.USE_PREF_SIZE
}
```

---

## RTL (Right-to-Left) Support

### CSS for RTL

```css
/* style-rtl.css */

.root:rtl {
  -fx-node-orientation: right-to-left;
}

/* Flip directional icons */
.button-back:rtl .glyph-icon {
  -fx-scale-x: -1;
}

/* Preserve LTR for scientific content */
.haplogroup-tree:rtl,
.chromosome-browser:rtl,
.sequence-id:rtl,
.snp-name:rtl {
  -fx-node-orientation: left-to-right;
}

/* Adjust split pane divider cursor */
.split-pane:rtl > .split-pane-divider {
  -fx-cursor: h-resize;
}
```

### ScalaFX RTL Setup

```scala
// In GenomeNavigatorApp.start()
private def applyLocaleOrientation(scene: Scene): Unit = {
  scene.root.value.nodeOrientation =
    if (I18n.isRTL) NodeOrientation.RightToLeft
    else NodeOrientation.LeftToRight
}

// Initial setup
applyLocaleOrientation(scene)

// Update on locale change
I18n.currentLocale.onChange { (_, _, _) =>
  applyLocaleOrientation(scene)

  // Reload RTL stylesheet if needed
  val stylesheets = scene.stylesheets
  stylesheets.removeIf(_.contains("style-rtl.css"))
  if (I18n.isRTL) {
    stylesheets.add(getClass.getResource("/style/style-rtl.css").toExternalForm)
  }
}
```

---

## Scientific Terms - Keep Universal

Some terms should **never** be translated:

### Universal Terms (haplogroups.properties)

```properties
# Haplogroup names - NEVER translate
# R1b-M269, H1a1, I1-M253, etc.

# Chromosome names
# Chr 1, Chr 2, X, Y, MT

# SNP names
# M269, P312, L21, etc.

# Gene names
# Keep as-is in all languages
```

### Universal Symbols

```properties
# Units - symbols are universal, labels may be localized
unit.centimorgans.symbol=cM
unit.coverage.symbol=x
unit.megabases.symbol=Mb
unit.basepairs.symbol=bp

# These remain the same in all locales
```

---

## Language Selector Implementation

```scala
// In SettingsDialog.scala

private def createLanguageSelector(): ComboBox[Locale] = {
  new ComboBox[Locale] {
    items = ObservableBuffer.from(I18n.supportedLocales)

    // Display language name in its own language
    cellFactory = _ => new ListCell[Locale] {
      item.onChange { (_, _, locale) =>
        text = Option(locale).map(_.getDisplayLanguage(locale)).getOrElse("")
      }
    }

    // Also set button cell
    buttonCell = new ListCell[Locale] {
      item.onChange { (_, _, locale) =>
        text = Option(locale).map(_.getDisplayLanguage(locale)).getOrElse("")
      }
    }

    // Initialize to current locale
    value = I18n.currentLocale.value

    // Handle selection change
    value.onChange { (_, _, newLocale) =>
      if (newLocale != null) {
        I18n.setLocale(newLocale)
        // Persist preference
        java.util.prefs.Preferences.userRoot()
          .node("com/decodingus")
          .put("locale", newLocale.toLanguageTag)
      }
    }
  }
}
```

---

## Testing Strategy

### Unit Tests

```scala
// src/test/scala/com/decodingus/i18n/I18nSpec.scala

class I18nSpec extends AnyFlatSpec with Matchers {

  "All message keys" should "exist in all supported locales" in {
    val baseBundle = ResourceBundle.getBundle("i18n.messages", Locale.ENGLISH)
    val baseKeys = baseBundle.getKeys.asScala.toSet

    val testLocales = Seq(Locale.GERMAN, new Locale("es"), Locale.FRENCH)

    testLocales.foreach { locale =>
      val bundle = ResourceBundle.getBundle("i18n.messages", locale)
      val localeKeys = bundle.getKeys.asScala.toSet

      val missing = baseKeys -- localeKeys
      withClue(s"Missing keys in ${locale.getDisplayLanguage}: ") {
        missing shouldBe empty
      }
    }
  }

  "Parameterized messages" should "format correctly in all locales" in {
    I18n.setLocale(Locale.ENGLISH)
    I18n.t("dashboard.subjects.count", 47) shouldBe "47 Subjects"

    I18n.setLocale(Locale.GERMAN)
    I18n.t("dashboard.subjects.count", 47) shouldBe "47 Probanden"
  }

  "Number formatting" should "use locale-specific separators" in {
    I18n.setLocale(Locale.ENGLISH)
    Formatters.formatNumber(1234567) shouldBe "1,234,567"

    I18n.setLocale(Locale.GERMAN)
    Formatters.formatNumber(1234567) shouldBe "1.234.567"
  }

  "RTL detection" should "identify RTL locales correctly" in {
    I18n.setLocale(Locale.ENGLISH)
    I18n.isRTL shouldBe false

    I18n.setLocale(new Locale("ar"))
    I18n.isRTL shouldBe true

    I18n.setLocale(new Locale("he"))
    I18n.isRTL shouldBe true
  }
}
```

### Visual Testing Checklist

- [ ] Test all screens in German (longest common translation)
- [ ] Verify no text truncation in buttons/labels
- [ ] Test RTL layout in Arabic
- [ ] Verify scientific terms remain LTR in RTL mode
- [ ] Check number/date formatting in each locale
- [ ] Verify locale persists across restarts

---

## Translation Workflow

### For Developers

1. Add new keys to `messages.properties` (English)
2. Run `sbt test` to catch missing translations
3. Coordinate with translators for new keys

### For Translators

1. Copy `messages.properties` to `messages_XX.properties`
2. Translate all values (keys remain in English)
3. Preserve placeholders like `{0}`, `{1}`
4. Test in application before submitting

### Key Guidelines for Translators

- Keep translations concise (UI space is limited)
- Preserve placeholder order unless grammar requires change
- Do not translate haplogroup names, SNP names, or scientific symbols
- Use formal/informal consistently (German: Sie vs du)
- Consider context (button vs heading vs description)

---

## Summary

| Aspect | Guideline |
|--------|-----------|
| Key naming | Hierarchical dot-notation (`nav.dashboard`) |
| Parameters | ICU MessageFormat with `{0}`, `{1}` |
| Reactive binding | Use `bind()` for dynamic locale updates |
| Layout | Flexible widths, no fixed text containers |
| RTL | Automatic via CSS, preserve LTR for science |
| Scientific terms | Keep universal (haplogroups, SNPs, units) |
| Testing | Automated key coverage + visual checks |
