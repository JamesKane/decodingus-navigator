ThisBuild / version := "0.1.0-SNAPSHOT"

ThisBuild / scalaVersion := "3.3.1"

// Require Java 17 LTS
ThisBuild / javacOptions ++= Seq("-source", "17", "-target", "17")

lazy val root = (project in file("."))
  .enablePlugins(JlinkPlugin)
  .settings(
    name := "DUNavigator",
    maintainer := "decodingus",

    // Main class for the application
    Compile / mainClass := Some("com.decodingus.ui.v2.NavigatorAppV2"),

    // Assembly settings for fat JAR
    assembly / assemblyJarName := s"DUNavigator-assembly-${version.value}.jar",
    assembly / mainClass := Some("com.decodingus.ui.v2.NavigatorAppV2"),
    // Enable Multi-Release JAR for Log4j2 compatibility with Java 9+
    assembly / packageOptions += Package.ManifestAttributes("Multi-Release" -> "true"),
    assembly / assemblyMergeStrategy := {
      case PathList("META-INF", "MANIFEST.MF") => MergeStrategy.discard
      case PathList("META-INF", "services", _*) => MergeStrategy.concat
      case PathList("META-INF", "native-image", _*) => MergeStrategy.discard
      case PathList("META-INF", "versions", _*) => MergeStrategy.first
      case PathList("META-INF", "LICENSE") => MergeStrategy.discard
      case PathList("META-INF", "LICENSE.txt") => MergeStrategy.discard
      case PathList("META-INF", "LICENSE.md") => MergeStrategy.discard
      case PathList("META-INF", "NOTICE") => MergeStrategy.discard
      case PathList("META-INF", "NOTICE.txt") => MergeStrategy.discard
      case PathList("META-INF", "NOTICE.md") => MergeStrategy.discard
      case PathList("META-INF", "DEPENDENCIES") => MergeStrategy.discard
      case PathList("META-INF", _*) => MergeStrategy.discard
      case "reference.conf" => MergeStrategy.concat
      case "application.conf" => MergeStrategy.concat
      case "module-info.class" => MergeStrategy.discard
      case x if x.endsWith(".proto") => MergeStrategy.first
      case x if x.endsWith(".class") => MergeStrategy.first
      case x if x.endsWith(".properties") => MergeStrategy.first
      case x if x.endsWith(".xml") => MergeStrategy.first
      case x if x.endsWith(".dtd") => MergeStrategy.first
      case x if x.endsWith(".xsd") => MergeStrategy.first
      case x if x.endsWith(".html") => MergeStrategy.first
      case x if x.endsWith(".txt") => MergeStrategy.first
      case x if x.contains("LICENSE") => MergeStrategy.discard
      case x if x.contains("NOTICE") => MergeStrategy.discard
      case x if x.contains("about.html") => MergeStrategy.discard
      case _ => MergeStrategy.first
    },

    resolvers += "Broad Institute" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release/",
    libraryDependencySchemes ++= Seq(
      "org.typelevel" %% "cats-kernel" % VersionScheme.Always,
      "org.typelevel" %% "cats-core" % VersionScheme.Always
    ),
    fork := true,

    // JVM options for forked process
    javaOptions ++= Seq(
      // Explicitly specify log4j config location to ensure ours is used
      "-Dlog4j.configurationFile=log4j2.xml",
      "-Dlog4j2.disable.jmx=true",
      // Open modules required by GATK/Spark reflection
      "--add-opens", "java.base/java.lang=ALL-UNNAMED",
      "--add-opens", "java.base/java.lang.invoke=ALL-UNNAMED",
      "--add-opens", "java.base/java.lang.reflect=ALL-UNNAMED",
      "--add-opens", "java.base/java.util=ALL-UNNAMED",
      "--add-opens", "java.base/java.io=ALL-UNNAMED",
      "--add-opens", "java.base/java.nio=ALL-UNNAMED",
      "--add-opens", "java.base/sun.nio.ch=ALL-UNNAMED"
    ),

    // Increase inline limit for Circe deriveCodec macros with deeply nested types
    scalacOptions += "-Xmax-inlines:256",
    libraryDependencies ++= {
      val jacksonVersion = "2.15.2"
      val javaFXVersion = "21.0.2"
      val osName = System.getProperty("os.name").toLowerCase match {
        case n if n.contains("mac") => "mac"
        case n if n.contains("win") => "win"
        case _ => "linux"
      }

      Seq(
        ("org.broadinstitute" % "gatk" % "4.6.2.0")
          .exclude("com.fasterxml.jackson.module", "jackson-module-scala_2.13")
          .exclude("com.fasterxml.jackson.core", "jackson-databind")
          .exclude("com.fasterxml.jackson.core", "jackson-core")
          .exclude("com.fasterxml.jackson.core", "jackson-annotations")
          .exclude("org.typelevel", "cats-kernel_2.13")
          .exclude("org.typelevel", "cats-core_2.13"),
        "com.fasterxml.jackson.module" %% "jackson-module-scala" % jacksonVersion,
        "com.fasterxml.jackson.core" % "jackson-databind" % jacksonVersion,
        "com.fasterxml.jackson.core" % "jackson-core" % jacksonVersion,
        "com.fasterxml.jackson.core" % "jackson-annotations" % jacksonVersion,
        "com.github.samtools" % "htsjdk" % "3.0.5",
        "co.fs2" %% "fs2-core" % "3.9.4",
        "io.circe" %% "circe-core" % "0.14.6",
        "io.circe" %% "circe-generic" % "0.14.6",
        "io.circe" %% "circe-parser" % "0.14.6",
        "org.scalafx" %% "scalafx" % "21.0.0-R32",
        "org.openjfx" % "javafx-base" % javaFXVersion classifier osName,
        "org.openjfx" % "javafx-controls" % javaFXVersion classifier osName,
        "org.openjfx" % "javafx-fxml" % javaFXVersion classifier osName,
        "org.openjfx" % "javafx-graphics" % javaFXVersion classifier osName,
        "org.openjfx" % "javafx-media" % javaFXVersion classifier osName,
        "org.openjfx" % "javafx-web" % javaFXVersion classifier osName,
        "com.typesafe" % "config" % "1.4.3",
        "com.softwaremill.sttp.client3" %% "core" % "3.9.7",
        "com.softwaremill.sttp.client3" %% "circe" % "3.9.7",
        // Database
        "com.h2database" % "h2" % "2.2.224",
        "com.zaxxer" % "HikariCP" % "5.1.0",
        "org.scalameta" %% "munit" % "1.0.0" % Test
      )
    },
    testFrameworks += new TestFramework("munit.Framework")
  )