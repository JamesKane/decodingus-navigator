package com.decodingus.analysis

import com.decodingus.workspace.model.ContigMetrics

case class CallableLociResult(
                               callableBases: Long,
                               contigAnalysis: List[ContigMetrics]
                             )
