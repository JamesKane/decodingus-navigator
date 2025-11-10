ThisBuild / version := "0.1.0-SNAPSHOT"

ThisBuild / scalaVersion := "3.3.1"

lazy val root = (project in file("."))
  .settings(
    name := "DUNavigator",
    resolvers += "Broad Institute" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release/",
    libraryDependencySchemes ++= Seq(
      "org.typelevel" %% "cats-kernel" % VersionScheme.Always,
      "org.typelevel" %% "cats-core" % VersionScheme.Always
    ),
    fork := true,
    libraryDependencies ++= {
      val jacksonVersion = "2.15.2"

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
        "org.scalafx" %% "scalafx" % "21.0.0-R32"
      )
    }
  )