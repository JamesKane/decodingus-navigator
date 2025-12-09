package com.decodingus.ancestry.model

import io.circe.Codec

/**
 * Defines a reference population for ancestry analysis.
 * Sub-continental granularity with ~25 populations from 1000 Genomes and HGDP/SGDP.
 */
case class Population(
  code: String,            // e.g., "CEU", "YRI", "CHB"
  name: String,            // "Northwestern European", "Yoruba", "Han Chinese"
  superPopulation: String, // "European", "African", "East Asian"
  region: String,          // Geographic region
  sampleCount: Int,        // Number of reference samples
  color: String            // Display color (hex)
) derives Codec.AsObject

object Population {

  // European populations
  val CEU = Population("CEU", "Northwestern European", "European", "Northern Europe", 99, "#0066CC")
  val FIN = Population("FIN", "Finnish", "European", "Northern Europe", 99, "#3399FF")
  val GBR = Population("GBR", "British", "European", "Northern Europe", 91, "#0044AA")
  val IBS = Population("IBS", "Iberian", "European", "Southern Europe", 107, "#66CCFF")
  val TSI = Population("TSI", "Tuscan", "European", "Southern Europe", 107, "#3366CC")

  // African populations
  val YRI = Population("YRI", "Yoruba", "African", "West Africa", 108, "#FF6600")
  val LWK = Population("LWK", "Luhya", "African", "East Africa", 99, "#FF9933")
  val ESN = Population("ESN", "Esan", "African", "West Africa", 99, "#CC6600")
  val MSL = Population("MSL", "Mende", "African", "West Africa", 85, "#FF8800")
  val GWD = Population("GWD", "Gambian", "African", "West Africa", 113, "#CC5500")

  // East Asian populations
  val CHB = Population("CHB", "Han Chinese", "East Asian", "East Asia", 103, "#00CC00")
  val JPT = Population("JPT", "Japanese", "East Asian", "East Asia", 104, "#33FF33")
  val KHV = Population("KHV", "Kinh Vietnamese", "East Asian", "Southeast Asia", 99, "#66FF66")
  val CHS = Population("CHS", "Southern Han Chinese", "East Asian", "East Asia", 105, "#00AA00")
  val CDX = Population("CDX", "Dai Chinese", "East Asian", "Southeast Asia", 93, "#009900")

  // South Asian populations
  val GIH = Population("GIH", "Gujarati", "South Asian", "South Asia", 103, "#9900CC")
  val PJL = Population("PJL", "Punjabi", "South Asian", "South Asia", 96, "#CC66FF")
  val BEB = Population("BEB", "Bengali", "South Asian", "South Asia", 86, "#9933FF")
  val STU = Population("STU", "Sri Lankan Tamil", "South Asian", "South Asia", 102, "#AA00DD")
  val ITU = Population("ITU", "Indian Telugu", "South Asian", "South Asia", 102, "#7700BB")

  // Americas populations (admixed)
  val MXL = Population("MXL", "Mexican", "Americas", "Central America", 64, "#CC0066")
  val PUR = Population("PUR", "Puerto Rican", "Americas", "Caribbean", 104, "#FF3399")
  val PEL = Population("PEL", "Peruvian", "Americas", "South America", 85, "#FF0066")
  val CLM = Population("CLM", "Colombian", "Americas", "South America", 94, "#DD0055")

  // HGDP populations for additional diversity
  val HGDP_Druze = Population("HGDP_Druze", "Druze", "West Asian", "Middle East", 47, "#996633")
  val HGDP_Palestinian = Population("HGDP_Palestinian", "Palestinian", "West Asian", "Middle East", 51, "#CC9966")
  val HGDP_Bedouin = Population("HGDP_Bedouin", "Bedouin", "West Asian", "Middle East", 48, "#AA8844")
  val HGDP_Papuan = Population("HGDP_Papuan", "Papuan", "Oceanian", "Oceania", 17, "#009999")
  val HGDP_Melanesian = Population("HGDP_Melanesian", "Melanesian", "Oceanian", "Oceania", 22, "#00BBBB")
  val HGDP_Yakut = Population("HGDP_Yakut", "Yakut", "Central Asian", "North Asia", 25, "#66CCCC")
  val HGDP_Maya = Population("HGDP_Maya", "Maya", "Native American", "Mesoamerica", 25, "#990033")
  val HGDP_Pima = Population("HGDP_Pima", "Pima", "Native American", "North America", 25, "#CC3366")
  val HGDP_Karitiana = Population("HGDP_Karitiana", "Karitiana", "Native American", "South America", 24, "#BB2255")

  /**
   * All reference populations used for ancestry estimation.
   */
  val All: List[Population] = List(
    // European
    CEU, FIN, GBR, IBS, TSI,
    // African
    YRI, LWK, ESN, MSL, GWD,
    // East Asian
    CHB, JPT, KHV, CHS, CDX,
    // South Asian
    GIH, PJL, BEB, STU, ITU,
    // Americas
    MXL, PUR, PEL, CLM,
    // West Asian (HGDP)
    HGDP_Druze, HGDP_Palestinian, HGDP_Bedouin,
    // Oceanian (HGDP)
    HGDP_Papuan, HGDP_Melanesian,
    // Central Asian (HGDP)
    HGDP_Yakut,
    // Native American (HGDP)
    HGDP_Maya, HGDP_Pima, HGDP_Karitiana
  )

  /**
   * Super-population groupings for summary display.
   */
  val SuperPopulations: Map[String, List[Population]] = All.groupBy(_.superPopulation)

  /**
   * Look up a population by its code.
   */
  def byCode(code: String): Option[Population] = All.find(_.code == code)

  /**
   * Get all populations for a super-population.
   */
  def bySuperPopulation(superPop: String): List[Population] =
    SuperPopulations.getOrElse(superPop, List.empty)

  /**
   * Total number of reference samples across all populations.
   */
  val totalSampleCount: Int = All.map(_.sampleCount).sum
}
