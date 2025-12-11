package com.decodingus.haplogroup.model

import io.circe.Codec

/**
 * Coordinate information for a variant in a specific reference build.
 */
case class VariantCoordinate(
  contig: String,
  position: Int,
  ref: String,
  alt: String
) derives Codec.AsObject

/**
 * Alias information for a variant, including common names, rsIds, and source attributions.
 */
case class VariantAliases(
  commonNames: List[String] = List.empty,
  rsIds: List[String] = List.empty,
  sources: Map[String, List[String]] = Map.empty
) derives Codec.AsObject

/**
 * Information about the haplogroup that this variant defines.
 */
case class DefiningHaplogroup(
  haplogroupId: Int,
  haplogroupName: String
) derives Codec.AsObject

/**
 * A named variant from the Decoding Us variant database.
 * This represents the public DTO returned by the variants API.
 *
 * @param variantId         Unique identifier for the variant
 * @param canonicalName     Primary/canonical name for the variant (e.g., "M269")
 * @param variantType       Type of variant (e.g., "SNP", "INDEL")
 * @param namingStatus      Current naming status
 * @param coordinates       Map of reference build to coordinate info (e.g., "GRCh38" -> coordinate)
 * @param aliases           Alternative names, rsIds, and source attributions
 * @param definingHaplogroup Optional haplogroup that this variant defines
 */
case class NamedVariant(
  variantId: Int,
  canonicalName: Option[String],
  variantType: String,
  namingStatus: String,
  coordinates: Map[String, VariantCoordinate],
  aliases: VariantAliases,
  definingHaplogroup: Option[DefiningHaplogroup] = None
) derives Codec.AsObject {

  /**
   * Get the coordinate for a specific reference build.
   */
  def coordinateFor(build: String): Option[VariantCoordinate] = coordinates.get(build)

  /**
   * Get all names for this variant (canonical + aliases).
   */
  def allNames: List[String] = {
    canonicalName.toList ++ aliases.commonNames
  }

  /**
   * Get a display name - canonical name if available, otherwise first common name.
   */
  def displayName: String = {
    canonicalName.orElse(aliases.commonNames.headOption).getOrElse(s"var_$variantId")
  }
}

/**
 * Metadata about the variant export file.
 */
case class VariantExportMetadata(
  generatedAt: String,
  variantCount: Int,
  fileSizeBytes: Long
) derives Codec.AsObject
