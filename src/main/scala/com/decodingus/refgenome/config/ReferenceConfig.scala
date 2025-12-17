package com.decodingus.refgenome.config

import io.circe.generic.semiauto.{deriveDecoder, deriveEncoder}
import io.circe.{Decoder, Encoder}

import java.nio.file.{Files, Path, Paths}

/**
 * Configuration for a single reference genome.
 *
 * @param build        The reference build name (e.g., "GRCh38", "GRCh37", "CHM13v2")
 * @param localPath    Optional user-specified local path to the reference FASTA (.fa.gz)
 * @param autoDownload Whether to automatically download if not found locally
 */
case class ReferenceGenomeConfig(
                                  build: String,
                                  localPath: Option[String] = None,
                                  autoDownload: Boolean = false
                                ) {
  /** Checks if the local path exists and is a valid file */
  def hasValidLocalPath: Boolean = localPath.exists { p =>
    val path = Paths.get(p)
    Files.exists(path) && Files.isRegularFile(path)
  }

  /** Gets the local path as a Path if valid */
  def getValidLocalPath: Option[Path] = localPath.flatMap { p =>
    val path = Paths.get(p)
    if (Files.exists(path) && Files.isRegularFile(path)) Some(path) else None
  }
}

object ReferenceGenomeConfig {
  implicit val encoder: Encoder[ReferenceGenomeConfig] = deriveEncoder
  implicit val decoder: Decoder[ReferenceGenomeConfig] = deriveDecoder
}

/**
 * Application-wide reference configuration.
 *
 * @param references           Map of build name to reference config
 * @param promptBeforeDownload Whether to prompt the user before downloading references
 * @param defaultCacheDir      Optional custom cache directory (defaults to ~/.decodingus/cache/references)
 */
case class ReferenceConfig(
                            references: Map[String, ReferenceGenomeConfig] = Map.empty,
                            promptBeforeDownload: Boolean = true,
                            defaultCacheDir: Option[String] = None
                          ) {
  /** Gets config for a specific build, or creates a default one */
  def getOrDefault(build: String): ReferenceGenomeConfig = {
    references.getOrElse(build, ReferenceGenomeConfig(build))
  }

  /** Updates config for a specific build */
  def withReference(config: ReferenceGenomeConfig): ReferenceConfig = {
    copy(references = references + (config.build -> config))
  }

  /** Gets the cache directory as a Path */
  def getCacheDir: Path = {
    defaultCacheDir.map(Paths.get(_)).getOrElse(
      Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "references")
    )
  }
}

object ReferenceConfig {
  implicit val encoder: Encoder[ReferenceConfig] = deriveEncoder
  implicit val decoder: Decoder[ReferenceConfig] = deriveDecoder

  /** Known reference builds with their download URLs */
  val knownBuilds: Map[String, String] = Map(
    "GRCh38" -> "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
    "GRCh37" -> "https://storage.googleapis.com/genomics-public-data/references/hg19/v0/Homo_sapiens_assembly19.fasta.gz",
    "CHM13v2" -> "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
  )

  /** Default configuration with all known builds */
  def default: ReferenceConfig = ReferenceConfig(
    references = knownBuilds.keys.map(b => b -> ReferenceGenomeConfig(b)).toMap,
    promptBeforeDownload = true
  )
}
