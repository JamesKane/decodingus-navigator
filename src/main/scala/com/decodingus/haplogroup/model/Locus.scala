package com.decodingus.haplogroup.model

import io.circe.Codec

enum LociType derives Codec.AsObject {
  case SNP, INDEL
}

case class LociCoordinate(
  position: Long,
  chromosome: String,
  ancestral: String,
  derived: String
) derives Codec.AsObject

case class Locus(
  name: String,
  loci_type: LociType,
  coordinates: Map[String, LociCoordinate]
) derives Codec.AsObject
