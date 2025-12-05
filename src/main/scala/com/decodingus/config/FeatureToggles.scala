package com.decodingus.config

import com.typesafe.config.ConfigFactory

object FeatureToggles {
  private val config = ConfigFactory.load("feature_toggles.conf")
  val pdsSubmissionEnabled: Boolean = config.getBoolean("pds-submission.enabled")
  val atProtocolEnabled: Boolean = config.hasPath("at-protocol.enabled") && config.getBoolean("at-protocol.enabled")

  object developerFeatures {
    private val devConfig = config.getConfig("developer-features")
    val saveJsonEnabled: Boolean = devConfig.getBoolean("save-json-enabled")
  }
}
