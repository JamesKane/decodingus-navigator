package com.decodingus.analysis

import com.decodingus.model.ContigSummary

case class CallableLociResult(
                               callableBases: Long,
                               contigAnalysis: List[ContigSummary]
                             )
