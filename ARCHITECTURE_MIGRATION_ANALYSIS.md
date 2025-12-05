## Architecture Analysis: Migration to Non-JVM Stack

The user is requesting alternatives to the current JVM-based architecture (Scala + GATK wrapper) due to persistent UI issues with JavaFX on macOS. They specifically asked about Rust.

### Current Architecture
*   **UI:** ScalaFX (JavaFX wrapper).
*   **Core Logic:**
    *   **WGS Metrics:** Wraps GATK `CollectWgsMetrics`.
    *   **Callable Loci:** Wraps GATK `CallableLoci`.
    *   **Haplogroup Calling:** Wraps GATK `HaplotypeCaller` (via `GatkHaplotypeCallerProcessor`) to call SNPs at specific tree loci, then scores them in Scala.
    *   **Liftover:** Wraps GATK `LiftoverVcf`.

### Feasibility of Migration

#### 1. WGS Metrics & Callable Loci (Coverage Analysis)
*   **Rust:** feasible. Tools like `mosdepth` (Nim, but fast) or libraries like `rust-htslib` / `noodles` can implemented this efficiently.
*   **Effort:** Moderate. Re-implementing the logic to match GATK's specific outputs (e.g., specific exclusion filters) requires care but is straightforward.

#### 2. Haplogroup Calling (The Blocker)
*   **Dependency:** The `HaplogroupProcessor.scala` relies on `GatkHaplotypeCallerProcessor`, which invokes GATK's `HaplotypeCaller`.
*   **Why GATK?** `HaplotypeCaller` performs local de-novo assembly to find variants. This is much more accurate than simple pileup-based calling, especially in repetitive regions or for indels.
*   **Rust Alternative:** `rust-htslib` allows pileup-based calling (reading bases at a position).
    *   *Pros:* Simple, fast.
    *   *Cons:* Less accurate than GATK HaplotypeCaller. Might miss variants or call false positives that GATK handles.
    *   *Mitigation:* For Haplogrouping, we are checking *known* SNPs. Simple pileup calling (checking if the read supports the derived allele) might be sufficient and is significantly easier to implement in Rust than a full HaplotypeCaller.

### Proposed "Hybrid" or "Pivot" Strategies

#### Option A: Rust CLI + Web UI (The "Rust" Request)
*   **Backend:** A new Rust CLI tool.
    *   Implements metrics using `rust-htslib`.
    *   Implements Haplogrouping using pileup-based calling (checking known loci) instead of full assembly.
*   **Frontend:** Tauri (Rust + React/HTML/JS) or just a local web server.
*   **Pros:** Single binary, small size, no JVM, modern UI.
*   **Cons:** Lower accuracy for haplogrouping (maybe acceptable?), rewrite of core logic.

#### Option B: Fix Current App (The "Path of Least Resistance")
*   **Problem:** JavaFX on macOS M1/M2 is finicky with modules.
*   **Solution:**
    1.  **Fat JAR:** Use `sbt assembly` to bundle everything. This often bypasses module path issues by putting everything on the classpath (though JavaFX 11+ warns against this, it usually works).
    2.  **JPackage:** Use `jpackage` to build a real `.dmg` / `.app` that bundles the correct JRE and JavaFX modules. This is the "correct" way to distribute JavaFX apps today.

### Recommendation for User
1.  **Short Term:** Try building a **Fat JAR** (`sbt assembly`). This is the quickest potential fix for the "run" issue without rewriting the app.
2.  **Long Term:** If GATK is a hard requirement, we can't easily leave the JVM without complicated containerization. If we can accept "Good Enough" pileup calling for haplogroups, a **Rust + Tauri** rewrite is extremely attractive for a consumer desktop app.

### Action Plan
1.  Add `sbt-assembly` plugin.
2.  Instruct user to run `sbt assembly`.
3.  Instruct user to run `java -jar target/.../DUNavigator.jar`.
