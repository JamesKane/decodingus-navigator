# Lexicon and Architecture Evolution Plan

## 1. Lexicon Updates

To capture the depth of sequencing statistics available in the `DUNavigator` (specifically from `WgsMetrics` and `LibraryStats`), we propose expanding the `sequenceData` record and introducing a dedicated `qualityMetrics` section.

```json
{
  "lexicon": 1,
  "id": "com.decodingus.atmosphere.workspace",
  "defs": {
    "main": {
      "type": "record",
      "description": "The root container for a Researcher's workspace, holding a pool of biosamples and defined projects.",
      "key": "tid",
      "record": {
        "type": "object",
        "required": ["samples", "projects"],
        "properties": {
          "samples": {
            "type": "array",
            "description": "The pool of all biosamples managed in this workspace.",
            "items": {
              "type": "ref",
              "ref": "#biosample"
            }
          },
          "projects": {
            "type": "array",
            "description": "Research projects grouping specific biosamples.",
            "items": {
              "type": "ref",
              "ref": "#project"
            }
          }
        }
      }
    },
    "biosample": {
      "type": "record",
      "description": "A record representing a biological sample and its associated sequencing metadata.",
      "key": "tid",
      "record": {
        "type": "object",
        "required": ["sampleAccession", "donorIdentifier", "centerName", "sequenceData"],
        "properties": {
          "sampleAccession": {
            "type": "string",
            "description": "Unique identifier for the sample (e.g., UUID from BGS)."
          },
          "donorIdentifier": {
            "type": "string",
            "description": "Identifier for the specimen donor within the user's context."
          },
          "description": {
            "type": "string",
            "description": "Human-readable description of the sample."
          },
          "centerName": {
            "type": "string",
            "description": "The name of the Sequencing Center or BGS Node."
          },
          "sex": {
            "type": "string",
            "description": "Biological sex of the donor.",
            "knownValues": ["Male", "Female", "Other", "Unknown"]
          },
          "sequenceData": {
            "type": "array",
            "description": "List of sequencing data entries, allowing for multiple alignments (e.g., GRCh38, chm13v2.0) or runs.",
            "items": {
              "type": "ref",
              "ref": "#sequenceData"
            }
          },
          "haplogroups": {
            "type": "ref",
            "ref": "#haplogroupAssignments",
            "description": "Y-DNA and mtDNA haplogroup assignments derived from the sequencing data."
          },
          "createdAt": {
            "type": "string",
            "format": "datetime"
          }
        }
      }
    },
    "project": {
      "type": "record",
      "description": "A genealogy or research project that aggregates multiple biosamples.",
      "key": "tid",
      "record": {
        "type": "object",
        "required": ["projectName", "administrator", "members"],
        "properties": {
          "projectName": {
            "type": "string",
            "description": "Name of the project (e.g., 'Smith Surname Project')."
          },
          "description": {
            "type": "string",
            "description": "Goals and scope of the research."
          },
          "administrator": {
            "type": "string",
            "description": "The DID or identifier of the researcher managing this project."
          },
          "members": {
            "type": "array",
            "description": "List of biosamples associated with this project.",
            "items": {
               "type": "string",
               "description": "Reference to a sampleAccession."
            }
          }
        }
      }
    },
    "haplogroupAssignments": {
      "type": "object",
      "description": "Container for paternal (Y-DNA) and maternal (mtDNA) haplogroup classifications.",
      "properties": {
        "yDna": {
          "type": "ref",
          "ref": "#haplogroupResult",
          "description": "The predicted Y-chromosome haplogroup (Paternal)."
        },
        "mtDna": {
          "type": "ref",
          "ref": "#haplogroupResult",
          "description": "The predicted Mitochondrial haplogroup (Maternal)."
        }
      }
    },
    "haplogroupResult": {
      "type": "object",
      "description": "Detailed scoring and classification result for a haplogroup.",
      "required": ["haplogroupName", "score"],
      "properties": {
        "haplogroupName": {
          "type": "string",
          "description": "The assigned haplogroup nomenclature (e.g., R-M269, H1a)."
        },
        "score": {
          "type": "float",
          "description": "Confidence score of the assignment."
        },
        "matchingSnps": {
          "type": "integer",
          "description": "Count of SNPs matching the defining mutations for this haplogroup."
        },
        "mismatchingSnps": {
          "type": "integer",
          "description": "Count of SNPs that contradict the assignment."
        },
        "ancestralMatches": {
          "type": "integer",
          "description": "Count of ancestral state matches."
        },
        "treeDepth": {
          "type": "integer",
          "description": "The depth of the assigned node in the phylogenetic tree."
        },
        "lineagePath": {
          "type": "array",
          "description": "The path from root to the assigned haplogroup (e.g., A -> ... -> R -> ... -> R-M269).",
          "items": {
            "type": "string"
          }
        }
      }
    },
    "sequenceData": {
      "type": "object",
      "description": "Raw sequencing run details and associated alignments.",
      "required": ["platformName", "testType", "files"],
      "properties": {
        "platformName": {
          "type": "string",
          "description": "Sequencing platform (e.g., ILLUMINA, PACBIO)."
        },
        "instrumentModel": {
          "type": "string",
          "description": "Specific instrument model (e.g., NovaSeq 6000)."
        },
        "testType": {
          "type": "string",
          "description": "Type of test (e.g., WGS, EXOME)."
        },
        "libraryLayout": {
          "type": "string",
          "description": "Paired-end or Single-end.",
          "knownValues": ["PAIRED", "SINGLE"]
        },
        "totalReads": {
          "type": "integer",
          "description": "Total number of reads."
        },
        "readLength": {
          "type": "integer",
          "description": "Average read length."
        },
        "meanInsertSize": {
          "type": "float",
          "description": "Mean insert size of the library."
        },
        "files": {
          "type": "array",
          "description": "Raw data files (e.g., FASTQs).",
          "items": {
            "type": "ref",
            "ref": "#fileInfo"
          }
        },
        "alignments": {
          "type": "array",
          "description": "List of alignments performed on this sequencing run.",
          "items": {
            "type": "ref",
            "ref": "#alignmentData"
          }
        }
      }
    },
    "alignmentData": {
      "type": "object",
      "description": "Details of a specific alignment (e.g., to GRCh38).",
      "required": ["referenceBuild", "aligner", "metrics"],
      "properties": {
        "referenceBuild": {
          "type": "string",
          "description": "Reference genome build (e.g., hg38, GRCh38)."
        },
        "aligner": {
          "type": "string",
          "description": "Tool used for alignment (e.g., BWA-MEM)."
        },
        "files": {
          "type": "array",
          "description": "Aligned data files (e.g., BAM, CRAM, VCF).",
          "items": {
            "type": "ref",
            "ref": "#fileInfo"
          }
        },
        "metrics": {
          "type": "ref",
          "ref": "#alignmentMetrics"
        }
      }
    },
    "alignmentMetrics": {
      "type": "object",
      "description": "Quality control metrics for the alignment.",
      "properties": {
        "genomeTerritory": {
          "type": "integer",
          "description": "The total number of bases in the reference genome territory."
        },
        "meanCoverage": {
          "type": "float",
          "description": "The mean coverage across the genome territory."
        },
        "medianCoverage": {
          "type": "float"
        },
        "sdCoverage": {
          "type": "float",
          "description": "Standard deviation of coverage."
        },
        "pctExcDupe": {
           "type": "float",
           "description": "Percentage of reads excluded due to duplication."
        },
        "pctExcMapq": {
           "type": "float",
           "description": "Percentage of reads excluded due to low mapping quality."
        },
        "pct10x": {
           "type": "float",
           "description": "Percentage of genome with at least 10x coverage."
        },
        "pct20x": {
           "type": "float",
           "description": "Percentage of genome with at least 20x coverage."
        },
        "pct30x": {
           "type": "float",
           "description": "Percentage of genome with at least 30x coverage."
        },
        "hetSnpSensitivity": {
           "type": "float",
           "description": "Sensitivity for detecting heterozygous SNPs."
        },
        "contigs": {
          "type": "array",
          "description": "Per-contig coverage statistics.",
          "items": {
             "type": "ref",
             "ref": "#contigMetrics"
          }
        }
      }
    },
    "contigMetrics": {
      "type": "object",
      "description": "Coverage analysis for a specific contig (chromosome).",
      "required": ["contigName", "callableBases"],
      "properties": {
        "contigName": {
          "type": "string",
          "description": "Name of the contig (e.g., chr1, 1)."
        },
        "callableBases": {
          "type": "integer",
          "description": "Number of bases deemed callable."
        },
        "meanCoverage": {
          "type": "float"
        },
        "poorMappingQuality": {
          "type": "integer",
          "description": "Number of bases with poor mapping quality."
        },
        "lowCoverage": {
           "type": "integer"
        },
        "noCoverage": {
           "type": "integer"
        }
      }
    },
    "fileInfo": {
      "type": "object",
      "description": "Metadata about a specific data file (FASTQ, BAM, etc.).",
      "required": ["fileName", "fileFormat", "location"],
      "properties": {
        "fileName": {
          "type": "string"
        },
        "fileSizeBytes": {
          "type": "integer"
        },
        "fileFormat": {
          "type": "string",
          "knownValues": ["FASTQ", "BAM", "CRAM", "VCF"]
        },
        "checksum": {
            "type": "string",
            "description": "SHA-256 or similar checksum."
        },
        "location": {
          "type": "string",
          "format": "uri",
          "description": "The URI where the file is stored (e.g., s3://..., ipfs://...)."
        }
      }
    }
  }
}
```

## 2. Architecture Evolution: Enabling the Federated Genealogy Network

The `DUNavigator` serves as the critical "Edge Ingestion Engine" for the Decoding Us ecosystem. Its primary role is to empower Citizen Scientist Researchers to process raw genomic data (BAM/CRAM) locally, extract standardized metadata (stats, haplogroups), and publish this "Genomic Passport" to a Personal Data Store (PDS).

### A. Data Hierarchy: Relational Project Management
We will adopt a relational model where `Biosamples` are independent entities that can be associated with one or more `Projects`.

1.  **Biosample (The Core Unit)**: The physical DNA source and its technical analysis (Genomic Passport). It exists independently of any project.
2.  **Project (The context)**: A "Virtual Folder" or grouping mechanism (e.g., "Smith Surname Project"). It contains metadata about the research goal and a list of references (`sampleAccession` IDs) to the Biosamples included in the study.
3.  **Subject (Participant)**: The donor entity. In this model, a Subject might be implicitly defined by the `donorIdentifier` on the Biosample or explicitly managed as a separate record if rich phenotype data is needed.

**Relationships:**
*   `Project` * -- * `Biosample` (A sample can be in the "Smith Surname" project AND the "R-U106" project).
*   `Biosample` 1 -- 1 `Analysis` (The passport is intrinsic to the sample processing event).

### B. Application Layer: The Researcher's Workbench

The UI/UX should evolve to support the workflow of a Project Administrator:

1.  **Project Manager Dashboard**:
    *   Create/Open local workspace/projects.
    *   Manage "Kits" (Subjects + Samples).
    *   **Feature**: "Bulk Import" from sequencing center manifests (e.g., FTDNA, Nebula, Dante).

2.  **Analysis & Classification**:
    *   **Haplogroup Caller**: Visualizing the placement of a sample on the Y-DNA/mtDNA tree.
    *   **QC Gatekeeper**: Automatically flagging samples that fail coverage thresholds (e.g., "Warning: < 10x coverage, Haplogroup prediction low confidence").

3.  **Local Workspace Persistence (Initial Stage)**:
    *   Instead of immediate PDS publication, the workbench will persist its entire state (all samples and projects) to a local `workspace.json` file.
    *   This file will be loaded on application startup and saved periodically or on explicit user action.
    *   This allows rapid prototyping and local management before full PDS integration.

### C. Enabling Federation & Matching (PDS Integration - Future)
While `DUNavigator` runs locally, it's designed to prepare the data for future integration with the broader network. The full PDS integration will involve:

*   **PDS Publication**: The Researcher "Publishes to PDS". This action signs the standardized JSON (Lexicon v1 - biosample records) and pushes it to the user's PDS.
*   **Discoverability**: Once in the PDS, the data becomes discoverable (permission-based) for the Federated Network.
*   **Standardization**: By adhering to the `com.decodingus.atmosphere.workspace` lexicon, all researchers speak the same language regarding coverage and haplogroups, enabling eventual federated queries and matching.
*   **Haplogroup Authority**: This tool acts as the initial "Authority" for haplogroup assignment before the data hits the network.
*   **IBD Readiness**: While IBD matching happens in the cloud/federation, this tool must ensure the *inputs* (e.g., phased VCFs or specific marker sets) are generated and validated during the processing step.
