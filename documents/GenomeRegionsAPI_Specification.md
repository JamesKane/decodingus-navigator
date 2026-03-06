# Genome Regions API Specification

**Status:** Pending Implementation
**Client Ready:** Yes (Navigator feature-flagged, disabled by default)
**Priority:** Medium

---

## Overview

The Navigator client requires a centralized API endpoint to provide curated genomic region metadata. This replaces the current approach of downloading from multiple external sources (ybrowse.org, GIAB, T2T S3) with a single authoritative source.

### Benefits

- **Single source of truth** - All region data from one curated endpoint
- **Simplified client** - No multi-source downloads, format parsing, or liftover logic
- **Better curation** - Human review before publishing; easy updates
- **CDN-ready** - Static JSON responses can be cached at edge
- **STR resolution** - CHM13v2.0 positions can be updated as they're verified

---

## Endpoint

```
GET /api/v1/genome-regions/{build}
```

### Path Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `build` | string | Yes | Reference genome build identifier |

### Supported Builds

| Build ID | Aliases | Description |
|----------|---------|-------------|
| `GRCh38` | `hg38` | GRCh38/hg38 (primary) |
| `GRCh37` | `hg19` | GRCh37/hg19 (legacy) |
| `CHM13v2` | `chm13`, `hs1`, `T2T-CHM13` | T2T CHM13v2.0 |

### Response Headers

```
Content-Type: application/json
Cache-Control: public, max-age=604800  # 7 days
ETag: "v2024.12.1"
```

---

## Response Schema

```json
{
  "build": "GRCh38",
  "version": "2024.12.1",
  "generatedAt": "2024-12-12T00:00:00Z",
  "chromosomes": {
    "chrY": {
      "length": 57227415,
      "centromere": { ... },
      "telomeres": { ... },
      "cytobands": [ ... ],
      "regions": { ... },
      "strMarkers": [ ... ]
    }
  }
}
```

### Top-Level Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `build` | string | Yes | Canonical build name (e.g., "GRCh38") |
| `version` | string | Yes | Semantic version for cache invalidation |
| `generatedAt` | ISO-8601 | Yes | Timestamp when data was generated |
| `chromosomes` | object | Yes | Map of chromosome name to region data |

### ChromosomeRegions Object

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `length` | integer | Yes | Total chromosome length in base pairs |
| `centromere` | Region | No | Centromere coordinates |
| `telomeres` | Telomeres | No | P-arm and Q-arm telomere regions |
| `cytobands` | Cytoband[] | Yes | Cytoband annotations for ideogram display |
| `regions` | YChromosomeRegions | No | Y-specific regions (only for chrY) |
| `strMarkers` | StrMarker[] | No | Named STR markers (only for chrY) |

### Region Object

```json
{
  "start": 10316944,
  "end": 10544039,
  "type": "Centromere",
  "modifier": 0.1
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `start` | integer | Yes | Start position (1-based, inclusive) |
| `end` | integer | Yes | End position (1-based, inclusive) |
| `type` | string | No | Region type classification |
| `modifier` | number | No | Quality modifier (1.0 = reliable, <1.0 = reduced confidence) |

### Telomeres Object

```json
{
  "p": { "start": 1, "end": 10000 },
  "q": { "start": 57217415, "end": 57227415 }
}
```

### Cytoband Object

```json
{
  "name": "q11.223",
  "start": 18500001,
  "end": 21000000,
  "stain": "gpos50"
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | Yes | Band name (e.g., "p36.33", "q11.21") |
| `start` | integer | Yes | Start position |
| `end` | integer | Yes | End position |
| `stain` | string | Yes | Giemsa stain pattern |

**Valid stain values:** `gneg`, `gpos25`, `gpos50`, `gpos75`, `gpos100`, `acen`, `gvar`, `stalk`

### YChromosomeRegions Object (chrY only)

```json
{
  "par1": { "start": 10001, "end": 2781479, "type": "PAR", "modifier": 0.5 },
  "par2": { "start": 56887903, "end": 57217415, "type": "PAR", "modifier": 0.5 },
  "xtr": [ ... ],
  "ampliconic": [ ... ],
  "palindromes": [ ... ],
  "heterochromatin": { ... },
  "xDegenerate": [ ... ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `par1` | Region | Pseudoautosomal region 1 (Yp) |
| `par2` | Region | Pseudoautosomal region 2 (Yq) |
| `xtr` | Region[] | X-transposed regions |
| `ampliconic` | Region[] | Ampliconic (high-copy) regions |
| `palindromes` | NamedRegion[] | Palindromic regions (P1-P8) |
| `heterochromatin` | Region | Yq12 heterochromatin |
| `xDegenerate` | Region[] | X-degenerate (stable single-copy) regions |

### NamedRegion Object

```json
{
  "name": "P1",
  "start": 24600000,
  "end": 25800000,
  "type": "Palindrome",
  "modifier": 0.4
}
```

### StrMarker Object

```json
{
  "name": "DYS389I",
  "start": 12684538,
  "end": 12684589,
  "period": 4,
  "verified": true,
  "note": null
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | Yes | Marker name (e.g., "DYS389I", "DYS456") |
| `start` | integer | Yes | Start position |
| `end` | integer | Yes | End position |
| `period` | integer | Yes | Repeat unit length in base pairs |
| `verified` | boolean | Yes | `true` if position manually verified for this build |
| `note` | string | No | Optional annotation (e.g., "Position estimated via liftover from GRCh38") |

---

## Quality Modifiers

The `modifier` field represents confidence in variant calls within each region type:

| Region Type | Modifier | Rationale |
|-------------|----------|-----------|
| X-degenerate | 1.0 | Stable, single-copy - gold standard |
| Normal | 1.0 | Standard callable region |
| PAR | 0.5 | Recombines with X chromosome |
| Palindrome | 0.4 | Gene conversion risk |
| XTR | 0.3 | 99% X-identical, contamination risk |
| Ampliconic | 0.3 | High copy number, mapping artifacts |
| Centromere | 0.1 | Nearly unmappable |
| Heterochromatin | 0.1 | Yq12 - essentially unmappable |

Modifiers combine multiplicatively for overlapping regions.

---

## Data Sources for Curation

The API backend should curate data from these authoritative sources:

### GRCh38

| Data Type | Source |
|-----------|--------|
| Cytobands | UCSC cytoBand table or ybrowse.org |
| Centromere | UCSC censat annotations |
| Telomeres | UCSC gap table |
| Y-PAR | GIAB genome-stratifications |
| Y-XTR | GIAB genome-stratifications |
| Y-Ampliconic | GIAB genome-stratifications |
| Y-Palindromes | ybrowse.org |
| Y-Heterochromatin | Manual curation (Yq12: 26673237-56887902) |
| Y-X-degenerate | Manual curation / T2T annotations |
| STR markers | ybrowse.org str_hg38.gff3 |

### GRCh37

Same sources with GRCh37/hg19 coordinates (native files available from ybrowse and GIAB).

### CHM13v2

| Data Type | Source |
|-----------|--------|
| Cytobands | T2T S3: chm13v2.0_cytobands_allchrs.bed |
| Centromere | T2T S3: chm13v2.0_censat_v2.1.bed |
| Y-PAR | GIAB genome-stratifications CHM13v2.0 |
| Y-XTR | GIAB genome-stratifications CHM13v2.0 |
| Y-Ampliconic | T2T S3: chm13v2.0Y_amplicons_v1.bed |
| Y-Palindromes | T2T S3: chm13v2.0Y_inverted_repeats_v1.bed |
| Y-Heterochromatin | T2T censat or manual (26637971-62122809) |
| Y-X-degenerate | T2T S3: chm13v2.0_chrXY_sequence_class_v1.bed |
| STR markers | Liftover from GRCh38 (mark `verified: false`) |

---

## Error Responses

### 404 Not Found

```json
{
  "error": "Unknown build",
  "message": "Build 'hg20' is not supported. Supported builds: GRCh38, GRCh37, CHM13v2",
  "supportedBuilds": ["GRCh38", "GRCh37", "CHM13v2"]
}
```

### 500 Internal Server Error

```json
{
  "error": "Internal error",
  "message": "Failed to load region data"
}
```

---

## Example Request/Response

### Request

```
GET /api/v1/genome-regions/GRCh38
Accept: application/json
```

### Response

```json
{
  "build": "GRCh38",
  "version": "2024.12.1",
  "generatedAt": "2024-12-12T00:00:00Z",
  "chromosomes": {
    "chrY": {
      "length": 57227415,
      "centromere": {
        "start": 10316944,
        "end": 10544039,
        "type": "Centromere",
        "modifier": 0.1
      },
      "telomeres": {
        "p": { "start": 1, "end": 10000 },
        "q": { "start": 57217415, "end": 57227415 }
      },
      "cytobands": [
        { "name": "p11.32", "start": 1, "end": 300000, "stain": "gneg" },
        { "name": "p11.31", "start": 300001, "end": 2800000, "stain": "gpos50" }
      ],
      "regions": {
        "par1": { "start": 10001, "end": 2781479, "type": "PAR", "modifier": 0.5 },
        "par2": { "start": 56887903, "end": 57217415, "type": "PAR", "modifier": 0.5 },
        "xtr": [
          { "start": 2781480, "end": 2963732, "type": "XTR", "modifier": 0.3 }
        ],
        "ampliconic": [
          { "start": 6171816, "end": 6681427, "type": "Ampliconic", "modifier": 0.3 }
        ],
        "palindromes": [
          { "name": "P1", "start": 24600000, "end": 25800000, "type": "Palindrome", "modifier": 0.4 }
        ],
        "heterochromatin": {
          "start": 26673237,
          "end": 56887902,
          "type": "Heterochromatin",
          "modifier": 0.1
        },
        "xDegenerate": [
          { "start": 2963733, "end": 6171815, "type": "XDegenerate", "modifier": 1.0 }
        ]
      },
      "strMarkers": [
        { "name": "DYS389I", "start": 12684538, "end": 12684589, "period": 4, "verified": true }
      ]
    }
  }
}
```

---

## Client Implementation Status

The Navigator client is ready to consume this API:

- **Feature toggle:** `genome-regions-api.enabled` (default: `false`)
- **Caching:** 7-day disk cache + in-memory cache
- **Fallback:** Bundled `grch38.json` for offline use
- **Legacy fallback:** File downloads from ybrowse/GIAB if API unavailable

### Enabling the API

Once deployed, enable in `feature_toggles.conf`:

```hocon
genome-regions-api {
  enabled = true
  base-url = "https://decoding-us.com/api/v1"
  cache-days = 7
  fallback-enabled = true
}
```

---

## Future Enhancements

1. **All chromosomes** - Currently focused on chrY; expand to include centromere/telomere data for all chromosomes
2. **MT-DNA regions** - Add mitochondrial region annotations
3. **Version diffing** - Endpoint to get changes between versions
4. **Bulk download** - Endpoint to download all builds in a single request

---

## Contact

For questions about this specification, contact the Navigator team.
