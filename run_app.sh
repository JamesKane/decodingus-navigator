#!/bin/bash

# Decoding-Us Navigator Launcher
# Runs the fat JAR created by 'sbt assembly' with necessary Java 21+ flags for JavaFX.

JAR_PATH="target/scala-3.3.1/DUNavigator-assembly-1.0.0-SNAPSHOT.jar"

if [ ! -f "$JAR_PATH" ]; then
    echo "Error: JAR file not found at $JAR_PATH"
    echo "Please run 'sbt assembly' first."
    exit 1
fi

echo "Starting Decoding-Us Navigator..."

java \
  --enable-native-access=ALL-UNNAMED \
  --add-modules=javafx.controls,javafx.fxml,javafx.graphics,javafx.media,javafx.web \
  --add-opens=javafx.base/com.sun.javafx.runtime=ALL-UNNAMED \
  --add-opens=javafx.controls/com.sun.javafx.scene.control.behavior=ALL-UNNAMED \
  --add-opens=javafx.controls/com.sun.javafx.scene.control=ALL-UNNAMED \
  --add-opens=javafx.base/com.sun.javafx.binding=ALL-UNNAMED \
  --add-opens=javafx.base/com.sun.javafx.event=ALL-UNNAMED \
  --add-opens=javafx.graphics/com.sun.javafx.stage=ALL-UNNAMED \
  --add-opens=javafx.graphics/com.sun.javafx.event=ALL-UNNAMED \
  --add-opens=javafx.graphics/com.sun.javafx.scene=ALL-UNNAMED \
  --add-opens=javafx.graphics/com.sun.javafx.sg.prism=ALL-UNNAMED \
  -jar "$JAR_PATH"
