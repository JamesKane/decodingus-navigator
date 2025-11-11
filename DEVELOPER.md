# DEVELOPER.md: Onboarding New Contributors

Welcome to the Decoding-Us Navigator project! This document aims to provide new contributors with an overview of the project's technology stack and identify areas where your contributions can make a significant impact.

## Project Overview

Decoding-Us Navigator is an Edge-computing companion to decoding.us.com. It enables local, privacy-preserving analysis of BAM/CRAM files for citizen scientists, leveraging GATK and providing a user-friendly interface. The application is built on the JVM for cross-platform compatibility.

## Technology Stack

### Language
*   **Scala 3:** The primary programming language, chosen for its conciseness, strong type system, and functional programming capabilities, running on the JVM.

### Build Tool
*   **SBT (Scala Build Tool):** Used for compiling, testing, and packaging the application.

### Core Libraries & Frameworks
*   **GATK (Genome Analysis Toolkit):** The core bioinformatics library used for various genomic analyses, such as CallableLoci.
*   **HTSJDK:** A Java API for accessing high-throughput sequencing data (BAM, CRAM, VCF, etc.). Used for reading and manipulating genomic files.
*   **FS2 (Functional Streams for Scala):** A functional streaming library used for asynchronous and efficient processing of large files and data streams.
*   **ScalaFX:** A Scala wrapper for JavaFX, used for building the graphical user interface (GUI) of the application.

### Configuration
*   **HOCON (Human-Optimized Config Object Notation):** Used for application configuration, including feature toggles (`feature_toggles.conf`).

### Data Formats
*   **JSON:** Used for data serialization, especially for summary statistics uploaded to the PDS.
*   **BED:** Output format for genomic regions from tools like GATK CallableLoci.
*   **BAM/CRAM:** Input formats for high-throughput sequencing data.

## Project Structure

The project follows a standard Scala/SBT directory structure:
*   `src/main/scala/com/decodingus/`: Contains the main Scala source code, organized by domain (e.g., `analysis`, `haplogroup`, `ui`).
*   `src/main/resources/`: Stores application resources, such as configuration files (`feature_toggles.conf`) and CSS (`style.css`).
*   `project/`: SBT build definitions.

## Areas for Contribution

We welcome contributions in various areas, from core bioinformatics logic to UI/UX improvements and testing.

### 1. Performance Improvements
*   **GATK Integration Optimization:** Explore ways to optimize GATK tool invocations, potentially through better parameter tuning or parallel execution strategies.
*   **Large File Processing:** Enhance the efficiency of reading and processing large BAM/CRAM files, especially for memory management and I/O operations.
*   **Parallelization:** Identify and implement further parallelization opportunities within the analysis pipelines to leverage multi-core processors more effectively.

### 2. Maintainability & Code Quality
*   **Code Refactoring:** Refactor existing code to improve clarity, reduce complexity, and adhere to functional programming principles where appropriate.
*   **Modularity:** Enhance the modularity of components to make them more independent and easier to test and maintain.
*   **Documentation:** Improve in-code documentation (Scaladoc) for complex functions, classes, and modules.

### 3. Test Coverage
*   **Unit Tests:** Expand unit test coverage for core logic components, especially in the `analysis` and `haplogroup` packages.
*   **Integration Tests:** Develop integration tests to ensure that different modules and external tool integrations (like GATK) work together seamlessly.
*   **UI Tests:** Implement tests for the ScalaFX user interface to ensure consistent behavior and responsiveness.

### 4. New Features & Enhancements
*   **Federation Integration:** Implement the "autosomal DNA matches with other researchers in the Federation" feature.
*   **Improved UI/UX:** Enhance the user interface with new visualizations, better feedback mechanisms, and more intuitive workflows.
*   **Additional Bioinformatics Tools:** Integrate other useful bioinformatics tools or analyses as identified by the community.
*   **Error Handling & Reporting:** Improve robust error handling and user-friendly error reporting.

## Getting Started

1.  **Clone the repository:** `git clone [repository-url]`
2.  **Install SBT:** Follow the instructions on the [SBT website](https://www.scala-sbt.org/download.html).
3.  **Open in IDE:** Import the project into your favorite Scala IDE (e.g., IntelliJ IDEA with Scala plugin).
4.  **Run the application:** `sbt run`
5.  **Run tests:** `sbt test`

We look forward to your contributions! If you have any questions, please don't hesitate to reach out.
