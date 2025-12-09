package com.decodingus.workspace.model

/**
 * A record representing a biological sample and its donor metadata.
 * This is a first-class record in the Atmosphere Lexicon (com.decodingus.atmosphere.biosample).
 *
 * Sequence and genotype data is referenced via AT URIs, not embedded.
 * This enables fine-grained CRUD operations without cascading changes.
 *
 * @param atUri                   The AT URI (at://did/collection/rkey) of this biosample record
 * @param meta                    Record metadata for tracking changes and sync
 * @param sampleAccession         Native identifier provided by the client for the biosample
 * @param donorIdentifier         Identifier for the specimen donor within the user's context
 * @param citizenDid              The Decentralized Identifier (DID) of the citizen/researcher who owns this record
 * @param description             Human-readable description of the sample
 * @param centerName              The name of the Sequencing Center or BGS Node
 * @param sex                     Biological sex of the donor (Male, Female, Other, Unknown)
 * @param haplogroups             Y-DNA and mtDNA haplogroup assignments derived from the sequencing data
 * @param sequenceRunRefs         AT URIs of sequence run records associated with this biosample
 * @param genotypeRefs            AT URIs of genotype data records (chip/array data) associated with this biosample
 * @param populationBreakdownRef  AT URI of the population/ancestry breakdown for this biosample
 * @param strProfileRef           AT URI of the Y-STR profile for this biosample (deprecated, use strProfileRefs)
 * @param strProfileRefs          AT URIs of Y-STR profile records for this biosample (supports multiple vendors)
 * @param ySnpPanelRefs           AT URIs of Y-DNA SNP panel result records (SNP packs, YSEQ panels, etc.)
 */
case class Biosample(
  atUri: Option[String],
  meta: RecordMeta,
  sampleAccession: String,
  donorIdentifier: String,
  citizenDid: Option[String] = None,
  description: Option[String] = None,
  centerName: Option[String] = None,
  sex: Option[String] = None,
  haplogroups: Option[HaplogroupAssignments] = None,
  sequenceRunRefs: List[String] = List.empty,
  genotypeRefs: List[String] = List.empty,
  populationBreakdownRef: Option[String] = None,
  strProfileRef: Option[String] = None,
  strProfileRefs: List[String] = List.empty,
  ySnpPanelRefs: List[String] = List.empty
)

object Biosample {
  /** Known sex values */
  val KnownSexValues: Set[String] = Set("Male", "Female", "Other", "Unknown")
}
