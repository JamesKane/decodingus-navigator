package com.decodingus.config

import com.typesafe.config.ConfigFactory

object FeatureToggles {
  private val config = ConfigFactory.load("feature_toggles.conf")
  val pdsSubmissionEnabled: Boolean = config.getBoolean("pds-submission.enabled")
}
