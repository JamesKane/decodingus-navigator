package com.decodingus.haplogroup.model

import io.circe.Codec

case class Locus(
                  name: String,
                  contig: String,
                  position: Long,
                  ref: String,
                  alt: String
                ) derives Codec.AsObject