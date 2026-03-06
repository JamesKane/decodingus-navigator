# Decoding-Us Navigator User Guide

Welcome to the **Decoding-Us Navigator**, your private, local companion for advanced genomic analysis. This application empowers you to explore your DNA data using professional-grade bioinformatics tools directly on your own computer, ensuring your genetic privacy while contributing to citizen science.

## Table of Contents
1. [Introduction](#introduction)
2. [System Requirements](#system-requirements)
3. [Installation & Setup](#installation--setup)
4. [Getting Started](#getting-started)
5. [Core Features](#core-features)
   - [Workspace Management](#workspace-management)
   - [Importing Sequencing Data](#importing-sequencing-data)
   - [Running Analyses](#running-analyses)
6. [Data Management & Privacy](#data-management--privacy)
7. [Advanced Usage](#advanced-usage)
8. [Troubleshooting](#troubleshooting)

---

## Introduction

Decoding-Us Navigator brings the power of the **Genome Analysis Toolkit (GATK)** to your desktop. Unlike cloud-based services where you must upload your raw DNA data, Navigator performs all heavy lifting locally. This "edge-computing" approach means:
- **Privacy First:** Your raw genomic files (BAM/CRAM) never leave your machine.
- **Data Sovereignty:** You own your data. Only optional, anonymized summaries are shared if you choose to connect to the Decoding-Us Federation.
- **Accessibility:** Complex command-line tools are wrapped in an easy-to-use graphical interface.

## System Requirements

To run Decoding-Us Navigator smoothly, your system should meet the following specifications:

- **Operating System:** macOS, Windows, or Linux.
- **Java Runtime:** **Java 17 LTS** or later is required.
- **Memory (RAM):** 
  - Minimum: 4 GB
  - Recommended: 8 GB or more (especially for large Whole Genome Sequencing files).
- **Disk Space:** Sufficient space for your BAM/CRAM files and downloaded reference genomes (approx. 5-10 GB for cached references).

## Installation & Setup

Currently, the application is distributed as a source project that can be built or run directly.

### Option 1: Running from Source (Recommended for testing)
If you have `sbt` (Scala Build Tool) installed:
1. Open a terminal in the project folder.
2. Run the application:
   ```bash
   sbt run
   ```

### Option 2: Building a Standalone App
You can create a standalone executable JAR file:
1. Run the build command:
   ```bash
   sbt assembly
   ```
2. The executable will be created in `target/scala-3.3.1/`. You can run it with Java:
   ```bash
   java -jar target/scala-3.3.1/DecodingUsNavigator-assembly-0.1.0-SNAPSHOT.jar
   ```

## Getting Started

### The Workbench
Upon launching, you will be greeted by the **Workbench**. This is your command center for organizing projects and research subjects.

- **Left Panel:** Navigation for Projects and Biosamples.
- **Center Panel:** Detailed views for the selected item (analysis results, charts, metrics).
- **Right Panel (Context):** Quick actions and details.

### First Launch
When you open Navigator for the first time, it will automatically initialize your local workspace database at `~/.decodingus/data/workspace.mv.db`. No manual configuration is needed.

## Core Features

### Workspace Management
Organize your research effectively:
- **Projects:** Create distinct projects to group related samples (e.g., "Family Study", "Ancient DNA").
- **Biosamples:** Create entries for individual research subjects.
- **Drag-and-Drop:** Easily move subjects between projects using drag-and-drop in the sidebar.

### Importing Sequencing Data
Navigator supports **BAM** and **CRAM** alignment files.
1. Select a Biosample in your workspace.
2. Click **"Import Sequence Run"** or drag-and-drop your file into the window.
3. The app automatically calculates a unique SHA-256 checksum to ensure data integrity and platform detection (Illumina, PacBio, Nanopore, etc.).

### Running Analyses
Once your data is imported, you can run various analysis modules. Results are cached, so re-running a module is instant if the data hasn't changed.

#### 1. Library Statistics
*Rapidly scan your file for quality metrics.*
- Checks read length distribution and insert sizes.
- Automatically detects the reference genome build (GRCh38, GRCh37, etc.).

#### 2. WGS Metrics
*Deep dive into coverage quality.*
- Visualizes coverage distribution (how much of the genome is covered and how deep).
- Calculates mean coverage and sensitivity thresholds (1x to 100x).

#### 3. Callable Loci
*Know exactly what you can trust.*
- Identifies "callable" bases—regions where the data is high quality enough to call variants reliably.
- Generates SVG visualizations of coverage gaps per chromosome.

#### 4. Haplogroup Determination
*Discover ancestral lineage.*
- **Y-DNA:** Traces paternal lineage using the Y chromosome.
- **MT-DNA:** Traces maternal lineage using mitochondrial DNA.
- **Private SNPs:** Identifies novel variants that define your unique branch on the tree.
- Supports multiple tree providers (FTDNA, DecodingUs).

#### 5. STR Profiles
*Short Tandem Repeat analysis.*
- Supports multi-panel STR profiling.
- Handles complex multi-allelic STRs.
- Compatible with various vendor formats.

## Data Management & Privacy

### Where is my data?
All application data is stored locally in your home directory under `~/.decodingus/`:
- **Database:** `~/.decodingus/data/` (Workspace structure, projects, samples)
- **Cache:** `~/.decodingus/cache/` (Analysis results, intermediate files)
- **References:** `~/.decodingus/cache/references/` (Downloaded genome builds)

### Reference Genomes
Navigator manages reference genomes automatically. It will download and cache standard builds (GRCh38, GRCh37, CHM13v2) as needed. You don't need to manually hunt for reference files.

### Cloud Integration (Experimental)
The application includes early support for the **AT Protocol** to sync workspaces.
*Note: This feature is currently in development and may be disabled by default or hidden behind feature flags.*
- **Goal:** Sync workspace metadata across devices via a Personal Data Store (PDS).
- **Privacy:** Even with cloud sync enabled, your raw BAM/CRAM files are **NEVER** uploaded. Only anonymized summary reports would be shared with your explicit consent.

## Advanced Usage

### Direct Database Access
Power users can query their workspace metadata directly using any H2-compatible database tool (like DBeaver).
- **JDBC URL:** `jdbc:h2:file:~/.decodingus/data/workspace`
- **Username:** `sa`
- **Password:** *(leave empty)*

**Note:** Ensure the Navigator app is closed before connecting externally to avoid file locks.

## Troubleshooting

**Q: The analysis is slow.**
A: WGS analysis is computationally intensive. Ensure you have assigned enough RAM to the application. If running via command line, you can increase memory: `sbt -mem 8192 run`.

**Q: I can't find my reference genome.**
A: Navigator downloads references on demand. If you are offline, ensure you have run an analysis at least once while online to cache the necessary files.

**Q: "Database is locked" error.**
A: This usually means another instance of Navigator is open, or an external tool is connected to the H2 database. Close other instances and try again.
