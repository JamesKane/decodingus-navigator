package com.decodingus.config

import io.circe.generic.semiauto.{deriveDecoder, deriveEncoder}
import io.circe.{Decoder, Encoder}

/**
 * Stored dimensions for a dialog window.
 */
case class DialogSize(width: Double, height: Double)

object DialogSize {
  implicit val decoder: Decoder[DialogSize] = deriveDecoder[DialogSize]
  implicit val encoder: Encoder[DialogSize] = deriveEncoder[DialogSize]
}

/**
 * User preferences for the application.
 * Stored in ~/.decodingus/config/user_preferences.json
 */
case class UserPreferences(
                            /** Y-DNA tree provider: "ftdna" or "decodingus" */
                            ydnaTreeProvider: String = "ftdna",

                            /** MT-DNA tree provider: "ftdna" or "decodingus" */
                            mtdnaTreeProvider: String = "ftdna",

                            /** UI locale as language tag (e.g., "en", "de", "es") */
                            locale: Option[String] = None,

                            /** UI theme: "dark" or "light" */
                            theme: Option[String] = None,

                            /** Saved dialog sizes keyed by dialog ID */
                            dialogSizes: Map[String, DialogSize] = Map.empty
                          )

object UserPreferences {
  val default: UserPreferences = UserPreferences()

  /** Valid tree provider values */
  val ValidTreeProviders: Set[String] = Set("ftdna", "decodingus")

  /** Display names for tree providers */
  def treeProviderDisplayName(provider: String): String = provider.toLowerCase match {
    case "ftdna" => "FTDNA (FamilyTreeDNA)"
    case "decodingus" => "Decoding-Us"
    case other => other
  }

  implicit val decoder: Decoder[UserPreferences] = deriveDecoder[UserPreferences]
  implicit val encoder: Encoder[UserPreferences] = deriveEncoder[UserPreferences]
}
