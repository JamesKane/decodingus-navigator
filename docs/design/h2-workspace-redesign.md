# H2 Workspace Database Redesign

## Executive Summary

This document describes a comprehensive redesign of the local Workspace caching layer, migrating from JSON-based persistence (`workspace.json`) to an embedded H2 database. The redesign optimizes for researchers managing large libraries of Subjects (100+ biosamples with associated sequence runs, alignments, and analysis artifacts).

### Goals

1. **Performance**: Sub-100ms queries for large workspaces (1000+ samples)
2. **Async PDS Sync**: Background synchronization with Personal Data Store
3. **File Cache Management**: Link tables for tracking analysis artifact state
4. **Collaboration-Ready**: Schema supports future cross-researcher deduplication (1K Genomes, HGDP)
5. **Transactional Integrity**: ACID guarantees for workspace operations

### Non-Goals

- Preserving alpha JSON cache (developers only)
- Implementing collaboration features (future work)
- Multi-user access (single desktop application)

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           DU Navigator Desktop                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────────────┐    ┌──────────────────┐    ┌──────────────────────┐   │
│  │ WorkbenchViewModel│◄──►│ WorkspaceService │◄──►│ H2WorkspaceRepository│   │
│  │ (Observable State)│    │ (Business Logic) │    │ (Data Access Layer)  │   │
│  └──────────────────┘    └────────┬─────────┘    └──────────┬───────────┘   │
│                                   │                         │               │
│                          ┌────────▼─────────┐    ┌──────────▼───────────┐   │
│                          │   SyncService    │    │  ~/.decodingus/      │   │
│                          │ (Async PDS Sync) │    │  workspace.h2.db     │   │
│                          └────────┬─────────┘    └──────────────────────┘   │
│                                   │                                         │
└───────────────────────────────────┼─────────────────────────────────────────┘
                                    │
                         ┌──────────▼──────────┐
                         │  Personal Data Store │
                         │  (AT Protocol)       │
                         └─────────────────────┘
```

### Component Responsibilities

| Component | Responsibility |
|-----------|---------------|
| `H2WorkspaceRepository` | CRUD operations, SQL queries, transaction management |
| `WorkspaceService` | Business logic, domain operations, cache coordination |
| `SyncService` | Async PDS push/pull, conflict detection, queue management |
| `WorkbenchViewModel` | Observable state for UI, delegates to services |

### Key Design Decisions

| Aspect | Decision |
|--------|----------|
| **Outgoing Sync** | Event-driven (immediate on edit) |
| **Incoming Sync** | Hourly poll (user can disable) |
| **Offline Support** | Indefinite - queue persists forever, no warnings |
| **Conflict UI** | Non-blocking status bar; user resolves when ready |
| **Canonical Lookup** | Schema prepared; implementation deferred to Q2 |

---

## Database Schema

### Core Entity Tables

```sql
-- Schema version tracking for migrations
CREATE TABLE schema_version (
    version INT PRIMARY KEY,
    applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    description VARCHAR(255)
);

-- Biosamples (Subjects)
CREATE TABLE biosample (
    id VARCHAR(36) PRIMARY KEY,           -- UUID
    at_uri VARCHAR(512),                  -- AT Protocol URI (nullable if not synced)
    sample_accession VARCHAR(255) NOT NULL,
    donor_identifier VARCHAR(255),
    description TEXT,
    center_name VARCHAR(255),
    biological_sex VARCHAR(20),           -- MALE, FEMALE, OTHER, UNKNOWN

    -- Haplogroup assignments (denormalized for performance)
    y_dna_haplogroup VARCHAR(255),
    y_dna_confidence DECIMAL(5,4),
    y_dna_tree_provider VARCHAR(50),
    y_dna_tree_version VARCHAR(50),
    mt_dna_haplogroup VARCHAR(255),
    mt_dna_confidence DECIMAL(5,4),
    mt_dna_tree_provider VARCHAR(50),
    mt_dna_tree_version VARCHAR(50),

    -- Record metadata
    version INT DEFAULT 1,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    last_modified_field VARCHAR(100),

    -- Sync state
    sync_status VARCHAR(20) DEFAULT 'NOT_SYNCED',  -- NOT_SYNCED, PENDING, SYNCED, MODIFIED, CONFLICT, ERROR
    at_cid VARCHAR(128),                  -- Content ID for version comparison
    last_synced_at TIMESTAMP,
    local_version INT DEFAULT 1,
    remote_version INT,

    -- Canonical sample tracking (for future collaboration)
    canonical_registry VARCHAR(50),        -- 1000GENOMES, HGDP, SGDP, ENA, NCBI
    canonical_accession VARCHAR(255),      -- Normalized accession

    CONSTRAINT uq_sample_accession UNIQUE (sample_accession)
);

CREATE INDEX idx_biosample_sync_status ON biosample(sync_status);
CREATE INDEX idx_biosample_canonical ON biosample(canonical_registry, canonical_accession);
CREATE INDEX idx_biosample_haplogroup_y ON biosample(y_dna_haplogroup);
CREATE INDEX idx_biosample_haplogroup_mt ON biosample(mt_dna_haplogroup);

-- Projects
CREATE TABLE project (
    id VARCHAR(36) PRIMARY KEY,
    at_uri VARCHAR(512),
    name VARCHAR(255) NOT NULL,
    description TEXT,
    administrator_did VARCHAR(255),

    -- Record metadata
    version INT DEFAULT 1,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    -- Sync state
    sync_status VARCHAR(20) DEFAULT 'NOT_SYNCED',
    at_cid VARCHAR(128),
    last_synced_at TIMESTAMP
);

-- Project-Biosample membership (many-to-many)
CREATE TABLE project_biosample (
    project_id VARCHAR(36) REFERENCES project(id) ON DELETE CASCADE,
    biosample_id VARCHAR(36) REFERENCES biosample(id) ON DELETE CASCADE,
    added_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (project_id, biosample_id)
);

-- Sequence Runs
CREATE TABLE sequence_run (
    id VARCHAR(36) PRIMARY KEY,
    at_uri VARCHAR(512),
    biosample_id VARCHAR(36) REFERENCES biosample(id) ON DELETE CASCADE,

    -- Platform details
    platform_name VARCHAR(50),             -- ILLUMINA, PACBIO, NANOPORE, etc.
    instrument_model VARCHAR(100),
    instrument_id VARCHAR(100),

    -- Test type
    test_type VARCHAR(50),                 -- WGS, WES, BIG_Y_500, MT_FULL_SEQUENCE, etc.

    -- Library info
    library_id VARCHAR(100),
    platform_unit VARCHAR(100),

    -- Sequencing metrics
    total_reads BIGINT,
    mapped_reads BIGINT,
    mean_insert_size INT,
    insert_size_std_dev INT,
    read_length INT,

    -- File reference (metadata only)
    source_filename VARCHAR(500),
    source_file_size BIGINT,
    source_file_checksum VARCHAR(128),
    source_file_format VARCHAR(20),        -- FASTQ, BAM, CRAM

    -- Capability flags
    supports_y_dna BOOLEAN DEFAULT FALSE,
    supports_mt_dna BOOLEAN DEFAULT FALSE,
    supports_ancestry BOOLEAN DEFAULT FALSE,

    -- Record metadata
    version INT DEFAULT 1,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    -- Sync state
    sync_status VARCHAR(20) DEFAULT 'NOT_SYNCED',
    at_cid VARCHAR(128),
    last_synced_at TIMESTAMP
);

CREATE INDEX idx_sequence_run_biosample ON sequence_run(biosample_id);
CREATE INDEX idx_sequence_run_platform ON sequence_run(platform_name);

-- Alignments
CREATE TABLE alignment (
    id VARCHAR(36) PRIMARY KEY,
    at_uri VARCHAR(512),
    sequence_run_id VARCHAR(36) REFERENCES sequence_run(id) ON DELETE CASCADE,

    -- Reference info
    reference_build VARCHAR(50),           -- GRCh38, GRCh37, T2T_CHM13, hg19, hg38
    aligner VARCHAR(100),
    aligner_version VARCHAR(50),
    variant_caller VARCHAR(100),
    variant_caller_version VARCHAR(50),

    -- Coverage metrics
    mean_coverage DECIMAL(10,2),
    median_coverage DECIMAL(10,2),
    coverage_std_dev DECIMAL(10,2),
    pct_bases_10x DECIMAL(5,2),
    pct_bases_20x DECIMAL(5,2),
    pct_bases_30x DECIMAL(5,2),

    -- Sex inference
    inferred_sex VARCHAR(20),
    sex_confidence DECIMAL(5,4),
    x_autosome_ratio DECIMAL(8,6),

    -- File reference (metadata only)
    alignment_filename VARCHAR(500),
    alignment_file_size BIGINT,
    alignment_file_checksum VARCHAR(128),
    alignment_file_format VARCHAR(20),     -- BAM, CRAM

    -- Record metadata
    version INT DEFAULT 1,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    -- Sync state
    sync_status VARCHAR(20) DEFAULT 'NOT_SYNCED',
    at_cid VARCHAR(128),
    last_synced_at TIMESTAMP
);

CREATE INDEX idx_alignment_sequence_run ON alignment(sequence_run_id);
CREATE INDEX idx_alignment_reference ON alignment(reference_build);

-- STR Profiles
CREATE TABLE str_profile (
    id VARCHAR(36) PRIMARY KEY,
    at_uri VARCHAR(512),
    biosample_id VARCHAR(36) REFERENCES biosample(id) ON DELETE CASCADE,

    -- Panel info
    panel_name VARCHAR(100),               -- Y-12, Y-37, Y-67, Y-111, Y-500, Y-700
    marker_count INT,
    provider VARCHAR(100),                 -- FTDNA, YSeq, YSTR-DB
    test_date DATE,

    -- Derivation info
    derivation_method VARCHAR(50),         -- WGS_DERIVED, VENDOR_PANEL, MANUAL_ENTRY
    source_run_id VARCHAR(36) REFERENCES sequence_run(id),

    -- Record metadata
    version INT DEFAULT 1,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    -- Sync state
    sync_status VARCHAR(20) DEFAULT 'NOT_SYNCED',
    at_cid VARCHAR(128),
    last_synced_at TIMESTAMP
);

CREATE INDEX idx_str_profile_biosample ON str_profile(biosample_id);

-- STR Marker Values (normalized)
CREATE TABLE str_marker_value (
    id VARCHAR(36) PRIMARY KEY,
    str_profile_id VARCHAR(36) REFERENCES str_profile(id) ON DELETE CASCADE,
    marker_name VARCHAR(50) NOT NULL,

    -- Value type discriminator
    value_type VARCHAR(20) NOT NULL,       -- SIMPLE, MULTI_COPY, COMPLEX

    -- Simple value
    simple_repeats INT,

    -- Multi-copy values (stored as JSON array for flexibility)
    multi_copy_values VARCHAR(100),        -- e.g., "[14,15]" for DYS385

    -- Complex value (alleles stored in separate table)
    raw_notation VARCHAR(100),

    -- Quality metrics
    quality_score DECIMAL(5,4),
    read_depth INT
);

CREATE INDEX idx_str_marker_profile ON str_marker_value(str_profile_id);
CREATE INDEX idx_str_marker_name ON str_marker_value(marker_name);

-- STR Complex Alleles (for palindromic markers)
CREATE TABLE str_allele (
    id VARCHAR(36) PRIMARY KEY,
    str_marker_value_id VARCHAR(36) REFERENCES str_marker_value(id) ON DELETE CASCADE,
    repeats INT NOT NULL,
    copy_count INT DEFAULT 1,
    designation VARCHAR(10)                -- 't', 'c', 'q' for tri/tetra/palindromic
);

-- Chip Profiles (SNP arrays)
CREATE TABLE chip_profile (
    id VARCHAR(36) PRIMARY KEY,
    at_uri VARCHAR(512),
    biosample_id VARCHAR(36) REFERENCES biosample(id) ON DELETE CASCADE,

    -- Vendor info
    vendor VARCHAR(50),                    -- 23ANDME, ANCESTRY, FTDNA, MYHERITAGE, LIVINGDNA
    test_type_code VARCHAR(20),
    chip_version VARCHAR(50),

    -- Quality metrics
    total_snps_called INT,
    no_call_rate DECIMAL(5,4),
    het_rate DECIMAL(5,4),

    -- DNA-specific counts
    y_marker_count INT,
    mt_marker_count INT,

    -- Import tracking
    import_date TIMESTAMP,
    source_file_hash VARCHAR(128),

    -- Record metadata
    version INT DEFAULT 1,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    -- Sync state
    sync_status VARCHAR(20) DEFAULT 'NOT_SYNCED',
    at_cid VARCHAR(128),
    last_synced_at TIMESTAMP
);

CREATE INDEX idx_chip_profile_biosample ON chip_profile(biosample_id);

-- Y-DNA SNP Panel Results
CREATE TABLE y_snp_panel_result (
    id VARCHAR(36) PRIMARY KEY,
    at_uri VARCHAR(512),
    biosample_id VARCHAR(36) REFERENCES biosample(id) ON DELETE CASCADE,

    -- Panel info
    panel_name VARCHAR(100),
    test_date DATE,
    vendor VARCHAR(100),

    -- Haplogroup result
    terminal_snp VARCHAR(100),
    haplogroup VARCHAR(255),
    confidence DECIMAL(5,4),

    -- Record metadata
    version INT DEFAULT 1,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    -- Sync state
    sync_status VARCHAR(20) DEFAULT 'NOT_SYNCED',
    at_cid VARCHAR(128),
    last_synced_at TIMESTAMP
);

-- Y-DNA SNP Calls (individual SNP results)
CREATE TABLE y_snp_call (
    id VARCHAR(36) PRIMARY KEY,
    panel_result_id VARCHAR(36) REFERENCES y_snp_panel_result(id) ON DELETE CASCADE,
    snp_name VARCHAR(100) NOT NULL,
    result VARCHAR(20) NOT NULL,           -- POSITIVE, NEGATIVE, NO_CALL
    grch38_position BIGINT,
    ancestral_allele CHAR(1),
    derived_allele CHAR(1)
);

CREATE INDEX idx_y_snp_call_panel ON y_snp_call(panel_result_id);
CREATE INDEX idx_y_snp_call_snp ON y_snp_call(snp_name);

-- Haplogroup Reconciliation (multi-run consensus)
CREATE TABLE haplogroup_reconciliation (
    id VARCHAR(36) PRIMARY KEY,
    at_uri VARCHAR(512),
    biosample_id VARCHAR(36) REFERENCES biosample(id) ON DELETE CASCADE,

    -- DNA type
    dna_type VARCHAR(10) NOT NULL,         -- Y_DNA, MT_DNA

    -- Consensus result
    consensus_haplogroup VARCHAR(255),
    consensus_confidence DECIMAL(5,4),
    compatibility_score DECIMAL(5,4),
    overall_compatibility VARCHAR(30),     -- COMPATIBLE, MINOR_DIVERGENCE, MAJOR_DIVERGENCE, INCOMPATIBLE

    -- Run tracking (JSON array of run IDs for flexibility)
    contributing_run_ids TEXT,             -- JSON array

    -- Record metadata
    version INT DEFAULT 1,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    -- Sync state
    sync_status VARCHAR(20) DEFAULT 'NOT_SYNCED',
    at_cid VARCHAR(128),
    last_synced_at TIMESTAMP
);

CREATE INDEX idx_reconciliation_biosample ON haplogroup_reconciliation(biosample_id);
CREATE INDEX idx_reconciliation_dna_type ON haplogroup_reconciliation(dna_type);

-- Run Haplogroup Calls (individual run contributions to reconciliation)
CREATE TABLE run_haplogroup_call (
    id VARCHAR(36) PRIMARY KEY,
    reconciliation_id VARCHAR(36) REFERENCES haplogroup_reconciliation(id) ON DELETE CASCADE,

    -- Source reference
    source_type VARCHAR(30) NOT NULL,      -- SEQUENCE_RUN, CHIP_PROFILE, Y_SNP_PANEL, STR_PROFILE
    source_id VARCHAR(36) NOT NULL,

    -- Call details
    haplogroup VARCHAR(255) NOT NULL,
    confidence DECIMAL(5,4),
    technology VARCHAR(30),                -- WGS, WES, BIG_Y, SNP_ARRAY, AMPLICON, STR_PANEL
    call_method VARCHAR(30),               -- SNP_PHYLOGENETIC, STR_PREDICTION, VENDOR_REPORTED
    supporting_snp_count INT,
    conflicting_snp_count INT,

    -- Weight for consensus calculation
    quality_tier INT                       -- 3=WGS, 2=BIG_Y, 1=Chip
);

CREATE INDEX idx_run_call_reconciliation ON run_haplogroup_call(reconciliation_id);
CREATE INDEX idx_run_call_source ON run_haplogroup_call(source_type, source_id);

-- SNP Conflicts within reconciliation
CREATE TABLE snp_conflict (
    id VARCHAR(36) PRIMARY KEY,
    reconciliation_id VARCHAR(36) REFERENCES haplogroup_reconciliation(id) ON DELETE CASCADE,

    snp_name VARCHAR(100) NOT NULL,
    resolution VARCHAR(30),                -- ACCEPT_MAJORITY, ACCEPT_HIGHER_QUALITY, ACCEPT_HIGHER_COVERAGE, UNRESOLVED, HETEROPLASMY
    resolved_value VARCHAR(20),
    notes TEXT
);

CREATE INDEX idx_snp_conflict_reconciliation ON snp_conflict(reconciliation_id);
```

### File Cache Link Tables

```sql
-- Analysis Artifact Cache Tracking
-- Links cached analysis outputs to their source alignments
CREATE TABLE analysis_artifact (
    id VARCHAR(36) PRIMARY KEY,
    alignment_id VARCHAR(36) REFERENCES alignment(id) ON DELETE CASCADE,

    -- Artifact type
    artifact_type VARCHAR(50) NOT NULL,    -- WGS_METRICS, CALLABLE_LOCI, HAPLOGROUP_VCF, WHOLE_GENOME_VCF, PRIVATE_VARIANTS

    -- Cache location
    cache_path VARCHAR(1000) NOT NULL,     -- Relative path within ~/.decodingus/cache/

    -- File metadata
    file_size BIGINT,
    file_checksum VARCHAR(128),
    file_format VARCHAR(20),

    -- Generation tracking
    generated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    generator_version VARCHAR(50),         -- GATK version, Navigator version
    generation_params TEXT,                -- JSON of parameters used

    -- Status
    status VARCHAR(20) DEFAULT 'AVAILABLE', -- AVAILABLE, IN_PROGRESS, STALE, DELETED
    stale_reason VARCHAR(255),             -- Why artifact is stale (e.g., reference changed)

    -- Dependencies (for invalidation)
    depends_on_source_checksum VARCHAR(128), -- Source BAM/CRAM checksum
    depends_on_reference_build VARCHAR(50)
);

CREATE INDEX idx_artifact_alignment ON analysis_artifact(alignment_id);
CREATE INDEX idx_artifact_type ON analysis_artifact(artifact_type);
CREATE INDEX idx_artifact_status ON analysis_artifact(status);

-- VCF Cache Metadata (detailed tracking for whole-genome VCFs)
CREATE TABLE vcf_cache (
    id VARCHAR(36) PRIMARY KEY,
    artifact_id VARCHAR(36) REFERENCES analysis_artifact(id) ON DELETE CASCADE,

    -- VCF-specific metadata
    reference_build VARCHAR(50) NOT NULL,
    variant_caller VARCHAR(100),
    variant_caller_version VARCHAR(50),
    gatk_version VARCHAR(50),

    -- Statistics
    total_variants INT,
    snp_count INT,
    indel_count INT,
    ti_tv_ratio DECIMAL(5,3),

    -- Contig coverage (JSON for flexibility)
    contig_stats TEXT,                     -- JSON: {"chrY": {"variants": 1234, "coverage": 30.5}, ...}

    -- Inferred metadata
    inferred_sex VARCHAR(20),
    sex_confidence DECIMAL(5,4),

    -- Index file tracking
    index_path VARCHAR(1000),
    index_checksum VARCHAR(128)
);

-- Callable Loci Cache (per-contig breakdown)
CREATE TABLE callable_loci_cache (
    id VARCHAR(36) PRIMARY KEY,
    artifact_id VARCHAR(36) REFERENCES analysis_artifact(id) ON DELETE CASCADE,

    -- Summary statistics
    total_callable_bases BIGINT,
    total_ref_n_bases BIGINT,
    total_no_coverage_bases BIGINT,
    total_low_coverage_bases BIGINT,

    -- Completion status
    is_complete BOOLEAN DEFAULT FALSE,
    completion_pct DECIMAL(5,2),

    -- Output files
    bed_file_path VARCHAR(1000),
    summary_table_path VARCHAR(1000),
    visualization_path VARCHAR(1000)       -- SVG
);

-- Callable Loci Contig Details
CREATE TABLE callable_loci_contig (
    id VARCHAR(36) PRIMARY KEY,
    callable_loci_cache_id VARCHAR(36) REFERENCES callable_loci_cache(id) ON DELETE CASCADE,

    contig_name VARCHAR(50) NOT NULL,
    callable_bases BIGINT,
    ref_n_bases BIGINT,
    no_coverage_bases BIGINT,
    low_coverage_bases BIGINT,
    callable_pct DECIMAL(5,2)
);

CREATE INDEX idx_callable_contig_cache ON callable_loci_contig(callable_loci_cache_id);

-- Source File Registry (tracking user's BAM/CRAM files)
CREATE TABLE source_file (
    id VARCHAR(36) PRIMARY KEY,
    alignment_id VARCHAR(36) REFERENCES alignment(id) ON DELETE SET NULL,

    -- File identity
    file_path VARCHAR(2000),               -- User's local path (may change)
    file_checksum VARCHAR(128) NOT NULL,   -- SHA-256 for stable identity
    file_size BIGINT,
    file_format VARCHAR(20),               -- BAM, CRAM

    -- Last known state
    last_verified_at TIMESTAMP,
    is_accessible BOOLEAN DEFAULT TRUE,

    -- Analysis state
    has_been_analyzed BOOLEAN DEFAULT FALSE,
    analysis_completed_at TIMESTAMP
);

CREATE INDEX idx_source_file_checksum ON source_file(file_checksum);
CREATE INDEX idx_source_file_alignment ON source_file(alignment_id);
```

### Sync Queue Tables

```sql
-- PDS Sync Queue (for async background sync)
-- Supports indefinite offline operation - no max_attempts limit
CREATE TABLE sync_queue (
    id VARCHAR(36) PRIMARY KEY,

    -- Target entity
    entity_type VARCHAR(50) NOT NULL,      -- BIOSAMPLE, PROJECT, SEQUENCE_RUN, ALIGNMENT, etc.
    entity_id VARCHAR(36) NOT NULL,

    -- Operation
    operation VARCHAR(20) NOT NULL,        -- CREATE, UPDATE, DELETE

    -- Queue state
    status VARCHAR(20) DEFAULT 'PENDING',  -- PENDING, IN_PROGRESS, COMPLETED, FAILED, RETRY
    priority INT DEFAULT 5,                -- 1-10, lower = higher priority

    -- Timing
    queued_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP,
    completed_at TIMESTAMP,

    -- Retry handling (no max_attempts - queue indefinitely for offline support)
    attempt_count INT DEFAULT 0,
    next_retry_at TIMESTAMP,
    last_error TEXT,

    -- Payload snapshot (JSON of entity state at queue time)
    payload_snapshot TEXT
);

CREATE INDEX idx_sync_queue_status ON sync_queue(status, priority, queued_at);
CREATE INDEX idx_sync_queue_entity ON sync_queue(entity_type, entity_id);

-- Sync History (audit trail)
CREATE TABLE sync_history (
    id VARCHAR(36) PRIMARY KEY,

    -- Target entity
    entity_type VARCHAR(50) NOT NULL,
    entity_id VARCHAR(36) NOT NULL,
    at_uri VARCHAR(512),

    -- Operation details
    operation VARCHAR(20) NOT NULL,
    direction VARCHAR(10) NOT NULL,        -- PUSH, PULL

    -- Result
    status VARCHAR(20) NOT NULL,           -- SUCCESS, FAILED, CONFLICT
    error_message TEXT,

    -- Timing
    started_at TIMESTAMP,
    completed_at TIMESTAMP,

    -- Version tracking
    local_version_before INT,
    local_version_after INT,
    remote_version_before INT,
    remote_version_after INT
);

CREATE INDEX idx_sync_history_entity ON sync_history(entity_type, entity_id);
CREATE INDEX idx_sync_history_time ON sync_history(completed_at);
```

### Future Collaboration Tables (Schema Only)

```sql
-- Canonical Sample Registry (for cross-researcher deduplication)
-- Populated when syncing samples matching known registries (1KG, HGDP, etc.)
CREATE TABLE canonical_sample (
    id VARCHAR(36) PRIMARY KEY,

    -- Registry identity
    registry VARCHAR(50) NOT NULL,         -- 1000GENOMES, HGDP, SGDP, ENA, NCBI
    canonical_accession VARCHAR(255) NOT NULL,

    -- Cross-references (populated from AppView)
    ena_accession VARCHAR(50),
    ncbi_biosample VARCHAR(50),

    -- Cached AppView merged values (refreshed on sync)
    appview_best_coverage DECIMAL(10,2),
    appview_y_haplogroup VARCHAR(255),
    appview_mt_haplogroup VARCHAR(255),
    appview_contributor_count INT,

    -- Local cache timestamp
    last_refreshed_at TIMESTAMP,

    CONSTRAINT uq_canonical UNIQUE (registry, canonical_accession)
);

CREATE INDEX idx_canonical_registry ON canonical_sample(registry, canonical_accession);

-- Link local biosamples to canonical samples
CREATE TABLE biosample_canonical_link (
    biosample_id VARCHAR(36) PRIMARY KEY REFERENCES biosample(id) ON DELETE CASCADE,
    canonical_sample_id VARCHAR(36) REFERENCES canonical_sample(id) ON DELETE SET NULL,

    -- Contribution metadata (what value our analysis adds)
    contributes_coverage BOOLEAN DEFAULT FALSE,
    contributes_y_haplogroup BOOLEAN DEFAULT FALSE,
    contributes_mt_haplogroup BOOLEAN DEFAULT FALSE,
    contributes_str_profile BOOLEAN DEFAULT FALSE,
    contributes_private_variants BOOLEAN DEFAULT FALSE,

    linked_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

---

## Repository Layer Design

### Base Repository Trait

```scala
package com.decodingus.workspace.repository

import java.sql.Connection
import scala.util.{Try, Using}

trait Repository[T, ID] {
  def findById(id: ID)(implicit conn: Connection): Option[T]
  def findAll()(implicit conn: Connection): List[T]
  def save(entity: T)(implicit conn: Connection): T
  def update(entity: T)(implicit conn: Connection): T
  def delete(id: ID)(implicit conn: Connection): Boolean
  def count()(implicit conn: Connection): Long
}

trait SyncableRepository[T, ID] extends Repository[T, ID] {
  def findPendingSync()(implicit conn: Connection): List[T]
  def findByAtUri(atUri: String)(implicit conn: Connection): Option[T]
  def updateSyncState(id: ID, state: SyncState)(implicit conn: Connection): Boolean
}
```

### Biosample Repository

```scala
package com.decodingus.workspace.repository

import com.decodingus.workspace.model.{Biosample, SyncState, SyncStatus}
import java.sql.{Connection, PreparedStatement, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

class BiosampleRepository extends SyncableRepository[Biosample, UUID] {

  override def findById(id: UUID)(implicit conn: Connection): Option[Biosample] = {
    val sql = """
      SELECT * FROM biosample WHERE id = ?
    """
    Using.resource(conn.prepareStatement(sql)) { stmt =>
      stmt.setString(1, id.toString)
      Using.resource(stmt.executeQuery()) { rs =>
        if (rs.next()) Some(mapRow(rs)) else None
      }
    }
  }

  def findBySampleAccession(accession: String)(implicit conn: Connection): Option[Biosample] = {
    val sql = """
      SELECT * FROM biosample WHERE sample_accession = ?
    """
    Using.resource(conn.prepareStatement(sql)) { stmt =>
      stmt.setString(1, accession)
      Using.resource(stmt.executeQuery()) { rs =>
        if (rs.next()) Some(mapRow(rs)) else None
      }
    }
  }

  def findByProject(projectId: UUID)(implicit conn: Connection): List[Biosample] = {
    val sql = """
      SELECT b.* FROM biosample b
      JOIN project_biosample pb ON b.id = pb.biosample_id
      WHERE pb.project_id = ?
      ORDER BY b.sample_accession
    """
    Using.resource(conn.prepareStatement(sql)) { stmt =>
      stmt.setString(1, projectId.toString)
      Using.resource(stmt.executeQuery()) { rs =>
        Iterator.continually(rs).takeWhile(_.next()).map(mapRow).toList
      }
    }
  }

  def findByHaplogroup(
    yDnaPattern: Option[String] = None,
    mtDnaPattern: Option[String] = None
  )(implicit conn: Connection): List[Biosample] = {
    val conditions = List(
      yDnaPattern.map(_ => "y_dna_haplogroup LIKE ?"),
      mtDnaPattern.map(_ => "mt_dna_haplogroup LIKE ?")
    ).flatten

    if (conditions.isEmpty) return List.empty

    val sql = s"""
      SELECT * FROM biosample
      WHERE ${conditions.mkString(" AND ")}
      ORDER BY sample_accession
    """
    Using.resource(conn.prepareStatement(sql)) { stmt =>
      var idx = 1
      yDnaPattern.foreach { p => stmt.setString(idx, p + "%"); idx += 1 }
      mtDnaPattern.foreach { p => stmt.setString(idx, p + "%"); idx += 1 }
      Using.resource(stmt.executeQuery()) { rs =>
        Iterator.continually(rs).takeWhile(_.next()).map(mapRow).toList
      }
    }
  }

  override def findPendingSync()(implicit conn: Connection): List[Biosample] = {
    val sql = """
      SELECT * FROM biosample
      WHERE sync_status IN ('NOT_SYNCED', 'MODIFIED', 'PENDING')
      ORDER BY updated_at
    """
    Using.resource(conn.prepareStatement(sql)) { stmt =>
      Using.resource(stmt.executeQuery()) { rs =>
        Iterator.continually(rs).takeWhile(_.next()).map(mapRow).toList
      }
    }
  }

  override def save(biosample: Biosample)(implicit conn: Connection): Biosample = {
    val sql = """
      INSERT INTO biosample (
        id, at_uri, sample_accession, donor_identifier, description,
        center_name, biological_sex, y_dna_haplogroup, y_dna_confidence,
        y_dna_tree_provider, y_dna_tree_version, mt_dna_haplogroup,
        mt_dna_confidence, mt_dna_tree_provider, mt_dna_tree_version,
        version, created_at, updated_at, sync_status, canonical_registry,
        canonical_accession
      ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    Using.resource(conn.prepareStatement(sql)) { stmt =>
      bindBiosample(stmt, biosample)
      stmt.executeUpdate()
    }
    biosample
  }

  // ... additional methods for update, delete, etc.

  private def mapRow(rs: ResultSet): Biosample = {
    // Map ResultSet to Biosample case class
    ???
  }

  private def bindBiosample(stmt: PreparedStatement, b: Biosample): Unit = {
    // Bind Biosample fields to PreparedStatement
    ???
  }
}
```

### Transaction Manager

```scala
package com.decodingus.workspace.repository

import java.sql.Connection
import javax.sql.DataSource
import scala.util.{Try, Using}

class TransactionManager(dataSource: DataSource) {

  def withTransaction[T](f: Connection => T): Try[T] = {
    Using(dataSource.getConnection()) { conn =>
      conn.setAutoCommit(false)
      try {
        val result = f(conn)
        conn.commit()
        result
      } catch {
        case e: Exception =>
          conn.rollback()
          throw e
      }
    }
  }

  def withReadOnly[T](f: Connection => T): Try[T] = {
    Using(dataSource.getConnection()) { conn =>
      conn.setReadOnly(true)
      f(conn)
    }
  }
}
```

---

## Async PDS Sync Architecture

### Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| **Outgoing Sync Trigger** | Event-driven (on edit) | Sync immediately when user makes changes |
| **Incoming Sync Frequency** | Hourly (optional) | Remote changes are infrequent; user can disable |
| **Offline Duration** | Indefinite | User may work offline permanently; no stale warnings |
| **Conflict Notification** | Status bar warning | Non-blocking; user resolves when ready |

### Sync Service Redesign

```scala
package com.decodingus.workspace.services

import com.decodingus.auth.User
import com.decodingus.workspace.repository._
import scala.concurrent.{ExecutionContext, Future}
import java.util.concurrent.{Executors, ScheduledExecutorService, TimeUnit}

class AsyncSyncService(
  transactionManager: TransactionManager,
  syncQueueRepository: SyncQueueRepository,
  pdsClient: PdsClient,
  conflictNotifier: ConflictNotifier
)(implicit ec: ExecutionContext) {

  private val scheduler: ScheduledExecutorService =
    Executors.newScheduledThreadPool(2)

  // Optional hourly incoming sync (can be disabled by user)
  private var incomingSyncEnabled: Boolean = true

  scheduler.scheduleWithFixedDelay(
    () => if (incomingSyncEnabled) pullRemoteChanges(),
    initialDelay = 60,
    delay = 3600,  // Hourly
    TimeUnit.SECONDS
  )

  def setIncomingSyncEnabled(enabled: Boolean): Unit = {
    incomingSyncEnabled = enabled
  }

  /**
   * Queue an entity for async sync to PDS.
   * Called immediately when user makes an edit.
   * Returns immediately; actual sync happens in background.
   */
  def queueForSync(
    entityType: String,
    entityId: String,
    operation: SyncOperation,
    priority: Int = 5
  ): Future[SyncQueueEntry] = Future {
    val entry = transactionManager.withTransaction { implicit conn =>
      syncQueueRepository.enqueue(SyncQueueEntry(
        entityType = entityType,
        entityId = entityId,
        operation = operation,
        priority = priority,
        status = SyncQueueStatus.Pending
      ))
    }.get

    // Trigger immediate processing for outgoing changes
    processOutgoingQueue()

    entry
  }

  /**
   * Process outgoing sync queue.
   * Called immediately after user edits (event-driven).
   * Queue persists indefinitely if offline - no timeout warnings.
   */
  private def processOutgoingQueue(): Unit = {
    transactionManager.withTransaction { implicit conn =>
      val pending = syncQueueRepository.findPendingBatch(batchSize = 10)
      pending.foreach { entry =>
        try {
          syncQueueRepository.markInProgress(entry.id)
          processSyncEntry(entry)
          syncQueueRepository.markCompleted(entry.id)
        } catch {
          case e: Exception =>
            handleSyncFailure(entry, e)
        }
      }
    }
  }

  /**
   * Pull remote changes from PDS.
   * Called hourly (if enabled) to detect conflicts.
   * User continues working; conflicts shown in status bar.
   */
  private def pullRemoteChanges(): Unit = {
    // Fetch remote state and detect conflicts
    // Conflicts are queued for user resolution, not blocking
    transactionManager.withTransaction { implicit conn =>
      val conflicts = detectRemoteConflicts()
      if (conflicts.nonEmpty) {
        conflictNotifier.notifyConflicts(conflicts)
      }
    }
  }

  private def processSyncEntry(entry: SyncQueueEntry)(implicit conn: Connection): Unit = {
    entry.operation match {
      case SyncOperation.Create => pushCreate(entry)
      case SyncOperation.Update => pushUpdate(entry)
      case SyncOperation.Delete => pushDelete(entry)
    }
  }

  /**
   * Handle sync failures with exponential backoff.
   * Queue persists indefinitely - user can work offline forever.
   */
  private def handleSyncFailure(entry: SyncQueueEntry, error: Exception)(
    implicit conn: Connection
  ): Unit = {
    // No max_attempts limit - queue indefinitely for offline support
    val backoff = Math.min(
      Math.pow(2, entry.attemptCount).toLong * 1000,
      3600000  // Cap at 1 hour between retries
    )
    syncQueueRepository.scheduleRetry(entry.id, backoff, error.getMessage)
  }

  def shutdown(): Unit = {
    scheduler.shutdown()
    scheduler.awaitTermination(30, TimeUnit.SECONDS)
  }
}
```

### Conflict Detection

```scala
package com.decodingus.workspace.services

sealed trait ConflictResolution
object ConflictResolution {
  case object KeepLocal extends ConflictResolution
  case object AcceptRemote extends ConflictResolution
  case class Merge(fields: Map[String, Any]) extends ConflictResolution
  case object RequireManual extends ConflictResolution
}

case class SyncConflict(
  entityType: String,
  entityId: String,
  localVersion: Int,
  remoteVersion: Int,
  localChanges: Set[String],       // Fields changed locally
  remoteChanges: Set[String],      // Fields changed remotely
  suggestedResolution: ConflictResolution
)

class ConflictDetector {

  // AppView-computed fields always take precedence
  private val appViewFields = Set(
    "y_dna_haplogroup", "y_dna_confidence",
    "mt_dna_haplogroup", "mt_dna_confidence"
  )

  def detectConflict(
    local: SyncableEntity,
    remote: SyncableEntity
  ): Option[SyncConflict] = {
    if (local.syncState.localVersion == remote.syncState.remoteVersion) {
      None // No conflict
    } else {
      val localChanges = detectChangedFields(local)
      val remoteChanges = detectChangedFields(remote)
      val overlapping = localChanges.intersect(remoteChanges)

      if (overlapping.isEmpty) {
        // Non-overlapping changes can be auto-merged
        Some(SyncConflict(
          entityType = local.entityType,
          entityId = local.id,
          localVersion = local.syncState.localVersion,
          remoteVersion = remote.syncState.remoteVersion.getOrElse(0),
          localChanges = localChanges,
          remoteChanges = remoteChanges,
          suggestedResolution = ConflictResolution.Merge(
            mergeFields(local, remote, localChanges, remoteChanges)
          )
        ))
      } else if (overlapping.subsetOf(appViewFields)) {
        // Only AppView fields conflict - accept remote
        Some(SyncConflict(
          entityType = local.entityType,
          entityId = local.id,
          localVersion = local.syncState.localVersion,
          remoteVersion = remote.syncState.remoteVersion.getOrElse(0),
          localChanges = localChanges,
          remoteChanges = remoteChanges,
          suggestedResolution = ConflictResolution.AcceptRemote
        ))
      } else {
        // True conflict - require manual resolution
        Some(SyncConflict(
          entityType = local.entityType,
          entityId = local.id,
          localVersion = local.syncState.localVersion,
          remoteVersion = remote.syncState.remoteVersion.getOrElse(0),
          localChanges = localChanges,
          remoteChanges = remoteChanges,
          suggestedResolution = ConflictResolution.RequireManual
        ))
      }
    }
  }
}
```

### Conflict UI: Status Bar Warning

Conflicts are presented via a non-blocking status bar warning. The user can continue working
and resolve conflicts when convenient.

```scala
package com.decodingus.workspace.services

import scalafx.beans.property.{ObjectProperty, BooleanProperty}
import scalafx.collections.ObservableBuffer

/**
 * Notifies the UI about sync conflicts via observable properties.
 * Non-blocking - user continues working until they choose to resolve.
 */
class ConflictNotifier {

  // Observable state for UI binding
  val hasConflicts: BooleanProperty = BooleanProperty(false)
  val conflictCount: ObjectProperty[Int] = ObjectProperty(0)
  val conflicts: ObservableBuffer[SyncConflict] = ObservableBuffer.empty

  def notifyConflicts(newConflicts: List[SyncConflict]): Unit = {
    Platform.runLater {
      conflicts.clear()
      conflicts ++= newConflicts
      conflictCount.value = newConflicts.size
      hasConflicts.value = newConflicts.nonEmpty
    }
  }

  def clearConflict(entityId: String): Unit = {
    Platform.runLater {
      conflicts.filterInPlace(_.entityId != entityId)
      conflictCount.value = conflicts.size
      hasConflicts.value = conflicts.nonEmpty
    }
  }
}
```

**Status Bar UI Design:**

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  Navigator                                                    [_] [□] [X]  │
├─────────────────────────────────────────────────────────────────────────────┤
│  [Menu Bar...]                                                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  [Main Content Area - User continues working normally]                      │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│  ⚠ 2 sync conflicts detected  [Review Conflicts]     │ ↻ 5 pending │ ✓ Online │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Conflict Resolution Dialog (opened when user clicks "Review Conflicts"):**

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  Sync Conflicts                                                    [X]      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  These items have been modified both locally and remotely.                  │
│  Your local changes are preserved until you choose how to resolve.          │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │ ⚠ VIK-002 (Biosample)                                                │  │
│  │   Local: Y-DNA haplogroup updated to R-Z284>BY3456                   │  │
│  │   Remote: Coverage updated to 45x by AppView                         │  │
│  │   Suggestion: Auto-merge (non-overlapping fields)                    │  │
│  │   [Accept Merge] [Keep Local] [Accept Remote] [View Details]         │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐  │
│  │ ⚠ ANC-007 (Biosample)                                                │  │
│  │   Local: mtDNA H1a                                                   │  │
│  │   Remote: mtDNA H1a1 (AppView haplogroup refinement)                 │  │
│  │   Suggestion: Accept remote (AppView refinement)                     │  │
│  │   [Accept Remote] [Keep Local] [View Details]                        │  │
│  └───────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
│  [Resolve All with Suggestions]                              [Close]        │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Key UX Principles:**

1. **Non-blocking**: User can ignore conflicts and continue working indefinitely
2. **Visible but unobtrusive**: Status bar shows conflict count without modal interruption
3. **Preserve local work**: Local changes are never lost without explicit user action
4. **Smart suggestions**: Auto-merge when safe; recommend AppView updates for haplogroups
5. **Batch resolution**: "Resolve All with Suggestions" for quick conflict clearing

---

## Database Connection Management

### H2 DataSource Configuration

```scala
package com.decodingus.workspace.db

import com.zaxxer.hikari.{HikariConfig, HikariDataSource}
import java.nio.file.{Files, Path, Paths}
import javax.sql.DataSource

object H2DataSource {

  private val dbDir: Path = Paths.get(
    System.getProperty("user.home"),
    ".decodingus"
  )

  private val dbPath: Path = dbDir.resolve("workspace")

  def createDataSource(): HikariDataSource = {
    // Ensure directory exists
    Files.createDirectories(dbDir)

    val config = new HikariConfig()
    config.setJdbcUrl(s"jdbc:h2:file:${dbPath};MODE=PostgreSQL;AUTO_SERVER=TRUE")
    config.setUsername("sa")
    config.setPassword("")
    config.setMaximumPoolSize(5)
    config.setMinimumIdle(1)
    config.setIdleTimeout(300000) // 5 minutes
    config.setConnectionTimeout(10000) // 10 seconds
    config.setMaxLifetime(1800000) // 30 minutes

    // H2-specific settings
    config.addDataSourceProperty("cachePrepStmts", "true")
    config.addDataSourceProperty("prepStmtCacheSize", "250")
    config.addDataSourceProperty("prepStmtCacheSqlLimit", "2048")

    new HikariDataSource(config)
  }
}
```

### Schema Migration

```scala
package com.decodingus.workspace.db

import java.sql.Connection
import scala.io.Source
import scala.util.Using

class SchemaMigrator(dataSource: javax.sql.DataSource) {

  private val migrations = List(
    1 -> "V1__initial_schema.sql",
    2 -> "V2__file_cache_tables.sql",
    3 -> "V3__sync_queue_tables.sql",
    4 -> "V4__canonical_sample_tables.sql"
  )

  def migrate(): Unit = {
    Using.resource(dataSource.getConnection()) { conn =>
      ensureSchemaVersionTable(conn)
      val currentVersion = getCurrentVersion(conn)

      migrations
        .filter(_._1 > currentVersion)
        .foreach { case (version, script) =>
          println(s"[DB] Applying migration $version: $script")
          applyMigration(conn, version, script)
        }
    }
  }

  private def ensureSchemaVersionTable(conn: Connection): Unit = {
    val sql = """
      CREATE TABLE IF NOT EXISTS schema_version (
        version INT PRIMARY KEY,
        applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        description VARCHAR(255)
      )
    """
    Using.resource(conn.createStatement()) { stmt =>
      stmt.execute(sql)
    }
  }

  private def getCurrentVersion(conn: Connection): Int = {
    val sql = "SELECT COALESCE(MAX(version), 0) FROM schema_version"
    Using.resource(conn.createStatement()) { stmt =>
      Using.resource(stmt.executeQuery(sql)) { rs =>
        rs.next()
        rs.getInt(1)
      }
    }
  }

  private def applyMigration(conn: Connection, version: Int, script: String): Unit = {
    val sql = Using.resource(
      Source.fromResource(s"db/migrations/$script")
    )(_.mkString)

    conn.setAutoCommit(false)
    try {
      Using.resource(conn.createStatement()) { stmt =>
        stmt.execute(sql)
      }

      Using.resource(conn.prepareStatement(
        "INSERT INTO schema_version (version, description) VALUES (?, ?)"
      )) { stmt =>
        stmt.setInt(1, version)
        stmt.setString(2, script)
        stmt.executeUpdate()
      }

      conn.commit()
    } catch {
      case e: Exception =>
        conn.rollback()
        throw new RuntimeException(s"Migration $version failed: ${e.getMessage}", e)
    } finally {
      conn.setAutoCommit(true)
    }
  }
}
```

---

## ViewModel Integration

### Updated WorkbenchViewModel

```scala
package com.decodingus.workspace

import com.decodingus.workspace.repository._
import com.decodingus.workspace.services._
import scalafx.beans.property.{ObjectProperty, BooleanProperty}
import scalafx.collections.ObservableBuffer
import scala.concurrent.ExecutionContext

class WorkbenchViewModel(
  transactionManager: TransactionManager,
  biosampleRepository: BiosampleRepository,
  projectRepository: ProjectRepository,
  syncService: AsyncSyncService
)(implicit ec: ExecutionContext) {

  // Observable state for UI binding
  val biosamples: ObservableBuffer[Biosample] = ObservableBuffer.empty
  val projects: ObservableBuffer[Project] = ObservableBuffer.empty
  val selectedBiosample: ObjectProperty[Option[Biosample]] = ObjectProperty(None)
  val syncInProgress: BooleanProperty = BooleanProperty(false)
  val pendingSyncCount: ObjectProperty[Int] = ObjectProperty(0)

  // Load initial state from H2
  def initialize(): Unit = {
    transactionManager.withReadOnly { implicit conn =>
      biosamples.clear()
      biosamples ++= biosampleRepository.findAll()

      projects.clear()
      projects ++= projectRepository.findAll()

      pendingSyncCount.value = syncQueueRepository.countPending()
    }
  }

  // Add biosample - triggers immediate sync attempt (event-driven)
  def addBiosample(biosample: Biosample): Unit = {
    transactionManager.withTransaction { implicit conn =>
      val saved = biosampleRepository.save(biosample)
      biosamples += saved

      // Queue and trigger immediate sync (event-driven on edit)
      // If offline, queued indefinitely until connectivity restored
      syncService.queueForSync(
        entityType = "BIOSAMPLE",
        entityId = saved.id.toString,
        operation = SyncOperation.Create
      )
    }
    updatePendingSyncCount()
  }

  // Update biosample - triggers immediate sync attempt (event-driven)
  def updateBiosample(biosample: Biosample): Unit = {
    transactionManager.withTransaction { implicit conn =>
      val updated = biosampleRepository.update(
        biosample.copy(
          meta = biosample.meta.updated("general"),
          syncState = biosample.syncState.copy(
            status = SyncStatus.Modified,
            localVersion = biosample.syncState.localVersion + 1
          )
        )
      )

      val idx = biosamples.indexWhere(_.id == updated.id)
      if (idx >= 0) biosamples.update(idx, updated)

      // Queue and trigger immediate sync (event-driven on edit)
      syncService.queueForSync(
        entityType = "BIOSAMPLE",
        entityId = updated.id.toString,
        operation = SyncOperation.Update
      )
    }
    updatePendingSyncCount()
  }

  // Force manual sync retry (user-triggered, e.g., after coming back online)
  def syncNow(): Unit = {
    syncInProgress.value = true
    syncService.processQueueNow().onComplete { _ =>
      Platform.runLater {
        syncInProgress.value = false
        updatePendingSyncCount()
        refreshFromDatabase()
      }
    }
  }

  private def updatePendingSyncCount(): Unit = {
    transactionManager.withReadOnly { implicit conn =>
      Platform.runLater {
        pendingSyncCount.value = syncQueueRepository.countPending()
      }
    }
  }

  private def refreshFromDatabase(): Unit = {
    transactionManager.withReadOnly { implicit conn =>
      val updatedBiosamples = biosampleRepository.findAll()
      Platform.runLater {
        biosamples.clear()
        biosamples ++= updatedBiosamples
      }
    }
  }
}
```

---

## Query Performance Optimizations

### Indexed Queries

```scala
// Fast haplogroup prefix search
def findByYDnaPrefix(prefix: String)(implicit conn: Connection): List[Biosample] = {
  val sql = """
    SELECT * FROM biosample
    WHERE y_dna_haplogroup LIKE ?
    ORDER BY y_dna_haplogroup
    LIMIT 1000
  """
  // Uses idx_biosample_haplogroup_y index
  ???
}

// Efficient project sample listing with pagination
def findByProjectPaginated(
  projectId: UUID,
  offset: Int,
  limit: Int
)(implicit conn: Connection): (List[Biosample], Int) = {
  val countSql = """
    SELECT COUNT(*) FROM project_biosample WHERE project_id = ?
  """
  val dataSql = """
    SELECT b.* FROM biosample b
    JOIN project_biosample pb ON b.id = pb.biosample_id
    WHERE pb.project_id = ?
    ORDER BY b.sample_accession
    LIMIT ? OFFSET ?
  """
  // Uses idx on project_biosample(project_id) and idx_biosample primary key
  ???
}

// Sync status dashboard query
def getSyncStatusSummary()(implicit conn: Connection): Map[SyncStatus, Int] = {
  val sql = """
    SELECT sync_status, COUNT(*) as cnt
    FROM biosample
    GROUP BY sync_status
  """
  // Uses idx_biosample_sync_status
  ???
}
```

### Batch Operations

```scala
// Bulk insert for CSV import
def insertBatch(biosamples: List[Biosample])(implicit conn: Connection): Int = {
  val sql = """
    INSERT INTO biosample (
      id, sample_accession, donor_identifier, description, biological_sex,
      sync_status, created_at, updated_at
    ) VALUES (?, ?, ?, ?, ?, 'NOT_SYNCED', CURRENT_TIMESTAMP, CURRENT_TIMESTAMP)
  """
  Using.resource(conn.prepareStatement(sql)) { stmt =>
    biosamples.foreach { b =>
      stmt.setString(1, b.id.toString)
      stmt.setString(2, b.sampleAccession)
      stmt.setString(3, b.donorIdentifier.orNull)
      stmt.setString(4, b.description.orNull)
      stmt.setString(5, b.biologicalSex.map(_.toString).orNull)
      stmt.addBatch()
    }
    stmt.executeBatch().sum
  }
}

// Bulk sync status update
def markAllAsPending(ids: List[UUID])(implicit conn: Connection): Int = {
  val sql = s"""
    UPDATE biosample
    SET sync_status = 'PENDING', updated_at = CURRENT_TIMESTAMP
    WHERE id IN (${ids.map(_ => "?").mkString(",")})
  """
  Using.resource(conn.prepareStatement(sql)) { stmt =>
    ids.zipWithIndex.foreach { case (id, idx) =>
      stmt.setString(idx + 1, id.toString)
    }
    stmt.executeUpdate()
  }
}
```

---

## Implementation Phases

### Phase 1: Core Database Layer ✅ COMPLETE
- [x] H2 DataSource setup and connection pooling (`Database.scala`)
- [x] Schema migration framework (`Migrator.scala`)
- [x] Core entity tables (`V001__initial_schema.sql`: biosample, project, project_member, sequence_run, alignment)
- [x] Basic repository implementations (`BiosampleRepository`, `ProjectRepository`, `SequenceRunRepository`, `AlignmentRepository`)
- [x] Transaction manager (`Transactor.scala`)
- [x] Comprehensive test coverage (127 repository tests)

### Phase 2: Full Entity Support ✅ COMPLETE
- [x] STR profile tables and repository (`V005__phase2_entity_tables.sql`, `StrProfileRepository`)
- [x] Chip profile tables and repository (`ChipProfileRepository`)
- [x] Y-DNA SNP panel tables and repository (`YSnpPanelRepository`)
- [x] Haplogroup reconciliation tables and repository (`HaplogroupReconciliationRepository`)
- [x] Project-biosample relationship management (`project_member` table, `ProjectRepository`)
- [x] Schema extensions for mixed panel imports:
  - `YSnpCall`: Added `endPosition`, `variantType` (SNP/INDEL), `orderedDate` for INDEL support
  - `StrMarkerValue`: Added `startPosition`, `endPosition`, `orderedDate` for genomic coordinate tracking
- [x] Backwards-compatible JSON codecs (reads old `position` field as `startPosition`)

### Phase 3: File Cache Integration ✅ COMPLETE
- [x] Analysis artifact tracking tables (`V002__file_cache_tables.sql`: `analysis_artifact`)
- [x] VCF cache metadata tables (`vcf_cache`)
- [x] Callable loci cache tables (`callable_loci_cache`, `callable_loci_contig`)
- [x] Source file registry (`source_file`, `SourceFileRepository`)
- [x] Coverage summary cache (`coverage_summary_cache`)
- [x] Unique constraint on artifact type per alignment (`V004__add_artifact_unique_constraint.sql`)
- [x] Cache invalidation logic (`CacheService` with `invalidateBySourceChecksum`, `invalidateByReferenceBuild`, `validateArtifact`, `validateAllArtifacts`, `verifyAllSourceFiles`, `cleanupMissingArtifacts`)

### Phase 4: Async PDS Sync ✅ COMPLETE
- [x] Sync queue tables (`V003__sync_queue_tables.sql`: `sync_queue`, `sync_history`, `sync_conflict`, `sync_settings`)
- [x] Sync queue repository (`SyncQueueRepository`)
- [x] Conflict detection and resolution (`SyncConflictRepository` with resolution actions)
- [x] Sync history and audit trail (`SyncHistoryRepository`)
- [x] Background sync worker service (`AsyncSyncService.scala` - queue processing, exponential backoff, scheduled polling)
- [x] UI sync status observable state (`ConflictNotifier.scala` - observable properties for UI binding)
- [x] Status bar UI component (`StatusBar.scala` - displays online/offline, pending count, conflicts; bound to ConflictNotifier)

### Phase 5: Future Collaboration Prep
- [ ] Canonical sample registry schema
- [ ] Biosample-canonical linking
- [ ] AppView query integration for contribution detection
- [ ] Deduplication UI hints

---

## Known Issues & Technical Notes

### H2 JSON Column Handling

**Issue**: H2's JSON data type maps to `byte[]` in Java, not `String`. Without proper handling:
- Writing strings directly creates JSON string literals (quoted values)
- Reading with `getString()` returns corrupted data

**Solution** (implemented in `Repository.scala`):
```scala
// Writing: Wrap JSON strings in JsonValue to convert to bytes
case class JsonValue(json: String)
def setParam(...) = value match {
  case JsonValue(j) => ps.setBytes(index, j.getBytes(UTF_8))
  ...
}

// Reading: Use getBytes() and decode as UTF-8
def getOptJsonString(rs: ResultSet, column: String): Option[String] =
  Option(rs.getBytes(column)).map(new String(_, UTF_8))
```

All repositories using JSON columns have been updated to use `JsonValue` wrapper.

### Pre-Phase 5 Checklist

Before starting Phase 5, address:
1. ~~**Background Sync Worker**: Service class implementing `AsyncSyncService` pattern from design doc~~ ✅ Complete
2. ~~**UI Sync Status**: Status bar integration showing pending/conflict counts~~ ✅ Complete (`StatusBar.scala` integrated into `GenomeNavigatorApp`)
3. ~~**Cache Invalidation**: Implement logic using `depends_on_source_checksum` and `depends_on_reference_build`~~ ✅ Complete
4. **Phase 2 Decision**: Determine if STR/Chip/SNP panel entities needed before collaboration features

---

## Testing Strategy

### Unit Tests ✅ IMPLEMENTED (198 tests passing)
- [x] Repository CRUD operations with in-memory H2 (all repositories)
- [x] Transaction rollback on error (`TransactorSpec`)
- [x] Conflict detection logic (`SyncConflictRepositorySpec`)
- [x] Migration script validation (`MigratorSpec`)
- [x] Database lifecycle (`DatabaseSpec`)
- [x] Sync queue operations (`SyncQueueRepositorySpec`)
- [x] Analysis artifact tracking (`AnalysisArtifactRepositorySpec`)
- [x] Source file registry (`SourceFileRepositorySpec`)
- [x] Async sync service (`AsyncSyncServiceSpec` - queue operations, stats, lifecycle)
- [x] Conflict notifier (`ConflictNotifierSpec` - observable state management)
- [x] Cache invalidation (`CacheServiceSpec` - artifact validation, source verification, cleanup)

### Integration Tests
- [ ] Full workflow: create → update → sync → conflict → resolve
- [ ] Large dataset performance (1000+ biosamples)
- [ ] Concurrent access simulation
- [ ] Cache invalidation scenarios

### Performance Benchmarks
- [ ] Query latency targets: <100ms for common operations
- [ ] Bulk import: 1000 samples in <10 seconds
- [ ] Sync queue throughput: 100 items/minute

---

## Migration Path

Since preserving the alpha JSON cache is out of scope:

1. On first launch with new version:
   - Initialize H2 database
   - Run schema migrations
   - Present empty workspace (fresh start)

2. Optional: Manual JSON import tool for developers
   - Read existing workspace.json
   - Parse and validate
   - Insert into H2 tables
   - Mark all as NOT_SYNCED

---

## Design Decisions

| Question | Decision | Implementation Notes |
|----------|----------|---------------------|
| **Sync Frequency** | Event-driven outgoing; hourly incoming (optional) | Outgoing: sync immediately on edit. Incoming: hourly poll for remote changes, user can disable entirely. |
| **Conflict UI** | Status bar warning, non-blocking | User continues working; conflicts shown in status bar with count. Resolution dialog opened on demand. |
| **Offline Duration** | Indefinite, no warnings | Queue persists forever. User may work offline permanently and never rejoin Federation. No stale data warnings. |
| **Canonical Lookup** | Deferred to Q2 next year | AppView canonical sample registry not in scope until Q2. Schema prepared but lookup logic deferred. |

---

## References

- [PDS Workbench Biosample Flow](../../decodingus/documents/proposals/pds-workbench-biosample-flow.md)
- [Atmosphere Lexicon Design](../../decodingus/documents/Atmosphere_Lexicon.md)
- [H2 Database Documentation](https://h2database.com/html/main.html)
- [HikariCP Configuration](https://github.com/brettwooldridge/HikariCP)
