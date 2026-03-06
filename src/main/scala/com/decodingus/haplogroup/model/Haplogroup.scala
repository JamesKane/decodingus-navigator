package com.decodingus.haplogroup.model

import io.circe.Codec

case class HaplogroupNode(
                           haplogroup_id: Long,
                           parent_id: Option[Long],
                           name: String,
                           is_root: Boolean,
                           loci: List[Locus],
                           children: List[Long]
                         ) derives Codec.AsObject

case class Haplogroup(
                       name: String,
                       parent: Option[String],
                       loci: List[Locus],
                       children: List[Haplogroup]
                     )

case class HaplogroupTree(
                           allNodes: Map[String, HaplogroupNode]
                         ) derives Codec.AsObject
