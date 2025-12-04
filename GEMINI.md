# Decoding-Us Navigator (DUNavigator)

## Project Overview
**Decoding-Us Navigator** is a desktop application designed for privacy-preserving, local genomic analysis. It serves as an edge-computing companion to [decoding.us.com](https://decoding.us.com), allowing citizen scientists and hobbyists to analyze BAM/CRAM files on their own machines without uploading raw data.

The application wraps complex command-line bioinformatics tools (specifically the Genome Analysis Toolkit - GATK) in a user-friendly GUI, enabling users to perform tasks like haplogroup classification, quality control (WGS metrics), and private SNP discovery.

## Technologies

*   **Language:** Scala 3.3.1
*   **Build Tool:** sbt (Scala Build Tool)
*   **UI Framework:** ScalaFX (Wrapper for JavaFX)
*   **Core Libraries:**
    *   **Bioinformatics:** GATK 4, HTSJDK
    *   **Streaming:** FS2 (Functional Streams for Scala)
    *   **JSON:** Circe, Jackson
    *   **Configuration:** Typesafe Config (HOCON)
    *   **HTTP:** sttp

## Build and Run

### Prerequisites
*   **Java Development Kit (JDK):** Required for Scala and JavaFX (Ensure compatibility with Scala 3.3.1).
*   **sbt:** The interactive build tool for Scala.

### Commands
*   **Run the Application:**
    ```bash
    sbt run
    ```
*   **Run Tests:**
    ```bash
    sbt test
    ```

## Project Structure

The project follows the standard sbt directory layout:

*   **`src/main/scala/com/decodingus/`**: Root package for source code.
    *   **`ui/`**: Contains the main application entry point (`GenomeNavigatorApp.scala`) and UI logic.
    *   **`analysis/`**: Core processors for bioinformatics tasks (e.g., `CallableLociProcessor`, `HaplogroupProcessor`).
    *   **`haplogroup/`**: Logic for Y-DNA and mtDNA haplogroup calling and tree data structures.
    *   **`model/`**: Domain models and data classes (e.g., `LibraryStats`, `WgsMetrics`).
    *   **`refgenome/`**: Manages reference genome downloads and caching.
    *   **`config/`**: Configuration loading (`FeatureToggles.scala`).
    *   **`pds/`**: Client for Personal Data Store interactions.
*   **`src/main/resources/`**: Non-code assets.
    *   `feature_toggles.conf`: Application configuration and feature flags.
    *   `style.css`: Stylesheet for the ScalaFX UI.

## Key Features & Components

1.  **Bioinformatics Pipeline:**
    *   **Library Stats:** Quickly estimates coverage, platform, and aligner from BAM/CRAM headers.
    *   **WGS Metrics:** Deep analysis of genome territory and coverage depth using GATK.
    *   **Callable Loci:** Identifies regions of the genome with sufficient coverage for reliable calling.
    *   **Haplogroup Calling:** Determines Y-DNA or mtDNA haplogroups using defined tree structures (DecodngUs or FTDNA).

2.  **UI (ScalaFX):**
    *   The interface is built using `StackPane`, `VBox`, and `HBox` layouts.
    *   It features a drag-and-drop area for input files and reactive progress bars for long-running analysis tasks.

3.  **Reference Management:**
    *   `ReferenceGateway` handles the resolution and downloading of reference genomes required for GATK analysis.

4.  **Concurrency:**
    *   Long-running tasks (like GATK processes) are offloaded from the JavaFX Application Thread using `javafx.concurrent.Task` and Scala futures/threads to keep the UI responsive.

## Development Conventions

*   **Functional Programming:** The project leverages functional patterns, particularly with the use of `FS2` for stream processing.
*   **Configuration:** Feature flags in `feature_toggles.conf` allow for enabling/disabling specific capabilities (e.g., PDS submission).
*   **Testing:** Unit and integration tests are expected for core logic. Refer to `DEVELOPER.md` for detailed contribution guidelines.
