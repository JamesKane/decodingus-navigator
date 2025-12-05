#!/bin/bash

# Script to package DUNavigator as a native macOS application (DMG)
# Prerequisite: JDK 14+ (with jpackage) installed and sbt.

APP_NAME="DUNavigator"
APP_VERSION="0.1.0"
MAIN_JAR="dunavigator.dunavigator-1.0.0-SNAPSHOT.jar"
MAIN_CLASS="com.decodingus.ui.GenomeNavigatorApp"
INPUT_DIR="target/universal/stage/lib"
OUTPUT_DIR="target/installer"

# JavaFX and Native Access Options
JAVA_OPTIONS=(
  "--enable-native-access=ALL-UNNAMED"
  "--add-modules=javafx.controls,javafx.fxml,javafx.graphics,javafx.media,javafx.web"
  "--add-opens=javafx.base/com.sun.javafx.runtime=ALL-UNNAMED"
  "--add-opens=javafx.controls/com.sun.javafx.scene.control.behavior=ALL-UNNAMED"
  "--add-opens=javafx.controls/com.sun.javafx.scene.control=ALL-UNNAMED"
  "--add-opens=javafx.base/com.sun.javafx.binding=ALL-UNNAMED"
  "--add-opens=javafx.base/com.sun.javafx.event=ALL-UNNAMED"
  "--add-opens=javafx.graphics/com.sun.javafx.stage=ALL-UNNAMED"
  "--add-opens=javafx.graphics/com.sun.javafx.event=ALL-UNNAMED"
  "--add-opens=javafx.graphics/com.sun.javafx.scene=ALL-UNNAMED"
  "--add-opens=javafx.graphics/com.sun.javafx.sg.prism=ALL-UNNAMED"
)

# Construct the java-options string for jpackage
# jpackage expects --java-options "opt1 opt2 ..." or multiple --java-options flags.
# We'll use multiple flags for safety.
JPACKAGE_OPTS=""
for opt in "${JAVA_OPTIONS[@]}"; do
  JPACKAGE_OPTS="$JPACKAGE_OPTS --java-options \"$opt\""
done

echo "Cleaning previous build..."
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

echo "Running sbt stage to prepare libraries..."
sbt stage

if [ ! -f "$INPUT_DIR/$MAIN_JAR" ]; then
    echo "Error: Main JAR not found at $INPUT_DIR/$MAIN_JAR"
    exit 1
fi

echo "Creating DMG installer with jpackage..."
# Note: We intentionally assume 'jpackage' is in the PATH.
# If using a specific JDK, set JAVA_HOME/bin/jpackage.

jpackage \
  --name "$APP_NAME" \
  --app-version "$APP_VERSION" \
  --input "$INPUT_DIR" \
  --main-jar "$MAIN_JAR" \
  --main-class "$MAIN_CLASS" \
  --type dmg \
  --dest "$OUTPUT_DIR" \
  --verbose \
  $JPACKAGE_OPTS

echo "Done. Installer should be in $OUTPUT_DIR"
