package com.decodingus.client

import com.decodingus.util.Logger
import io.circe.*
import io.circe.generic.auto.*
import sttp.client3.*
import sttp.client3.circe.*

import scala.concurrent.{ExecutionContext, Future}

/**
 * Metadata resolved from the European Nucleotide Archive for a biological sample.
 *
 * @param sampleAccession ENA/INSDC accession (SAMN/SAMEA)
 * @param sampleAlias     Human-readable alias (e.g., HG02759)
 * @param description     Free-text description from the archive
 * @param sex             Biological sex if available
 * @param centerName      Submitting center/institution
 * @param taxId           NCBI taxonomy ID (9606 for Homo sapiens)
 * @param population      Population code (e.g., PEL) — typically from 1KG samples
 * @param populationName  Full population name (e.g., Peruvian in Lima, Peru)
 * @param superPopulation Continental superpopulation code (e.g., AMR)
 */
case class EnaSampleMetadata(
  sampleAccession: String,
  sampleAlias: Option[String] = None,
  description: Option[String] = None,
  sex: Option[String] = None,
  centerName: Option[String] = None,
  taxId: Option[String] = None,
  population: Option[String] = None,
  populationName: Option[String] = None,
  superPopulation: Option[String] = None
)

/**
 * Metadata for an ENA study/project.
 * Designed for future Project sync support.
 */
case class EnaStudyMetadata(
  studyAccession: String,
  studyTitle: Option[String] = None,
  studyDescription: Option[String] = None,
  centerName: Option[String] = None
)

/**
 * Client for the European Nucleotide Archive (ENA) Portal API and the
 * EMBL-EBI BioSamples API. Used to resolve sample metadata from public archives
 * given a sample alias (e.g., HG02759) or accession (e.g., SAMN00009587).
 *
 * Resolution strategy:
 *   1. ENA Portal API — fast search by alias or accession, returns basic fields
 *   2. BioSamples API — rich attribute lookup (sex, population) by accession
 *
 * Designed for extension: `resolveStudySamples` supports future bulk Project sync from ENA.
 */
object EnaClient {

  private val log = Logger[EnaClient.type]
  private val backend = HttpClientFutureBackend()

  private val EnaPortalBase = "https://www.ebi.ac.uk/ena/portal/api"
  private val BioSamplesBase = "https://www.ebi.ac.uk/biosamples/samples"

  // ENA Portal API fields to request
  private val PortalFields = "accession,sample_alias,description,center_name,tax_id"

  // ── ENA Portal API DTOs ──────────────────────────────────────────────

  private case class EnaPortalRow(
    accession: Option[String] = None,
    sample_alias: Option[String] = None,
    description: Option[String] = None,
    center_name: Option[String] = None,
    tax_id: Option[String] = None
  )

  // ── BioSamples API DTOs ──────────────────────────────────────────────

  private case class BioSampleCharacteristicValue(text: Option[String] = None)
  private case class BioSamplesResponse(
    accession: Option[String] = None,
    name: Option[String] = None,
    characteristics: Option[Map[String, List[BioSampleCharacteristicValue]]] = None
  )

  // ── Public API ───────────────────────────────────────────────────────

  /**
   * Resolve metadata for a single sample by alias or accession.
   * Tries ENA Portal first for basic fields, then enriches from BioSamples API.
   */
  def resolveSample(identifier: String)(implicit ec: ExecutionContext): Future[Option[EnaSampleMetadata]] = {
    log.info(s"Resolving sample metadata for: $identifier")

    // Try as alias first, then as accession
    resolveViaPortal(identifier).flatMap {
      case Some(base) =>
        // Enrich with BioSamples attributes (sex, population)
        enrichFromBioSamples(base)
      case None =>
        // identifier might already be an accession — try BioSamples directly
        fetchBioSample(identifier).map {
          case Some(enriched) => Some(enriched)
          case None =>
            log.info(s"No ENA/BioSamples record found for: $identifier")
            None
        }
    }.recover { case e: Exception =>
      log.error(s"Failed to resolve sample $identifier: ${e.getMessage}")
      None
    }
  }

  /**
   * Resolve metadata for all samples in a study/project.
   * Designed for future Project sync — fetches all sample accessions in one call,
   * then enriches each from BioSamples.
   */
  def resolveStudySamples(studyAccession: String)(implicit ec: ExecutionContext): Future[List[EnaSampleMetadata]] = {
    log.info(s"Resolving all samples for study: $studyAccession")

    val query = s"""study_accession="$studyAccession""""
    val request = basicRequest
      .get(uri"$EnaPortalBase/search?result=sample&query=$query&fields=$PortalFields&format=json&limit=0")
      .response(asJson[List[EnaPortalRow]])
      .readTimeout(scala.concurrent.duration.Duration(30, "s"))

    request.send(backend).flatMap { response =>
      response.body match {
        case Right(rows) =>
          log.info(s"Found ${rows.size} samples in study $studyAccession")
          val bases = rows.map(portalRowToMetadata)
          // Enrich in parallel (bounded)
          Future.traverse(bases)(enrichFromBioSamples).map(_.flatten)
        case Left(error) =>
          log.error(s"ENA Portal search failed for study $studyAccession: $error")
          Future.successful(List.empty)
      }
    }.recover { case e: Exception =>
      log.error(s"Failed to resolve study samples: ${e.getMessage}")
      List.empty
    }
  }

  /**
   * Resolve study-level metadata.
   * Designed for future Project sync.
   */
  def resolveStudy(studyAccession: String)(implicit ec: ExecutionContext): Future[Option[EnaStudyMetadata]] = {
    val query = s"""study_accession="$studyAccession""""
    val fields = "study_accession,study_title,study_description,center_name"
    val request = basicRequest
      .get(uri"$EnaPortalBase/search?result=study&query=$query&fields=$fields&format=json&limit=1")
      .response(asJson[List[Map[String, String]]])

    request.send(backend).map { response =>
      response.body.toOption.flatMap(_.headOption).map { row =>
        EnaStudyMetadata(
          studyAccession = row.getOrElse("study_accession", studyAccession),
          studyTitle = row.get("study_title").filter(_.nonEmpty),
          studyDescription = row.get("study_description").filter(_.nonEmpty),
          centerName = row.get("center_name").filter(_.nonEmpty)
        )
      }
    }.recover { case _ => None }
  }

  // ── Internal ─────────────────────────────────────────────────────────

  /**
   * Search ENA Portal API by alias, falling back to accession search.
   */
  private def resolveViaPortal(identifier: String)(implicit ec: ExecutionContext): Future[Option[EnaSampleMetadata]] = {
    // Try alias first
    searchPortal(s"""sample_alias="$identifier"""").flatMap {
      case Some(meta) => Future.successful(Some(meta))
      case None =>
        // Try as accession
        searchPortal(s"""accession="$identifier"""")
    }
  }

  private def searchPortal(query: String)(implicit ec: ExecutionContext): Future[Option[EnaSampleMetadata]] = {
    val request = basicRequest
      .get(uri"$EnaPortalBase/search?result=sample&query=$query&fields=$PortalFields&format=json&limit=1")
      .response(asJson[List[EnaPortalRow]])

    request.send(backend).map { response =>
      response.body.toOption.flatMap(_.headOption).map(portalRowToMetadata)
    }.recover { case _ => None }
  }

  private def portalRowToMetadata(row: EnaPortalRow): EnaSampleMetadata =
    EnaSampleMetadata(
      sampleAccession = row.accession.getOrElse(""),
      sampleAlias = row.sample_alias.filter(_.nonEmpty),
      description = row.description.filter(_.nonEmpty),
      centerName = row.center_name.filter(_.nonEmpty),
      taxId = row.tax_id.filter(_.nonEmpty)
    )

  /**
   * Enrich base metadata with BioSamples API attributes (sex, population, etc.).
   */
  private def enrichFromBioSamples(base: EnaSampleMetadata)(implicit ec: ExecutionContext): Future[Option[EnaSampleMetadata]] = {
    if (base.sampleAccession.isEmpty) return Future.successful(Some(base))

    fetchBioSample(base.sampleAccession).map {
      case Some(enriched) =>
        Some(base.copy(
          sampleAlias = base.sampleAlias.orElse(enriched.sampleAlias),
          sex = enriched.sex.orElse(base.sex),
          population = enriched.population,
          populationName = enriched.populationName,
          superPopulation = enriched.superPopulation,
          description = base.description.orElse(enriched.description)
        ))
      case None =>
        Some(base)
    }
  }

  /**
   * Fetch rich sample attributes from the BioSamples API.
   */
  private def fetchBioSample(accession: String)(implicit ec: ExecutionContext): Future[Option[EnaSampleMetadata]] = {
    val request = basicRequest
      .get(uri"$BioSamplesBase/$accession")
      .response(asJson[BioSamplesResponse])

    request.send(backend).map { response =>
      response.body.toOption.map { bs =>
        val chars = bs.characteristics.getOrElse(Map.empty)
        EnaSampleMetadata(
          sampleAccession = bs.accession.getOrElse(accession),
          sampleAlias = bs.name.filter(_.nonEmpty),
          sex = extractCharacteristic(chars, "sex", "Sex", "biological sex"),
          population = extractCharacteristic(chars, "population", "Population"),
          populationName = extractCharacteristic(chars, "population description", "Population Description",
            "population name", "Population Name"),
          superPopulation = extractCharacteristic(chars, "superpopulation", "Superpopulation",
            "superpopulation code", "super population")
        )
      }
    }.recover { case _ => None }
  }

  /**
   * Extract a characteristic value from BioSamples, trying multiple key variants.
   * BioSamples attribute names are inconsistent across submissions.
   */
  private def extractCharacteristic(
    chars: Map[String, List[BioSampleCharacteristicValue]],
    keys: String*
  ): Option[String] = {
    keys.view
      .flatMap(key => chars.get(key))
      .flatMap(_.headOption)
      .flatMap(_.text)
      .filter(_.nonEmpty)
      .headOption
  }
}
