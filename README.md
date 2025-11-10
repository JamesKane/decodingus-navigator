# **Decoding-Us Navigator: Callable Loci Analyzer (Scala 3 / SBT)**

This project serves as the foundation for the Decoding-us.com Edge computing service, designed to perform local analysis of BAM/CRAM files and securely contribute anonymized summary statistics to the user's Personal Data Store (PDS).  
The core application replaces the initial Ruby scripts, using the native JVM ecosystem (Scala/Java) for robust performance and deployment alongside the GATK toolkit.

## **🚀 1\. Project Setup (SBT)**

This project uses Scala 3 and is managed by SBT (Scala Build Tool).

### **1.1 Dependencies (build.sbt)**

Key dependencies required for this project include:

| Dependency | Purpose |
| :---- | :---- |
| org.scala-lang | Scala 3 Standard Library |
| org.broadinstitute.gatk | GATK Core Libraries (for CallableLoci functionality) |
| com.github.samtools | HTSJDK (for BAM/CRAM reading and index/dict generation) |
| com.outr / fs2 | Asynchronous file I/O and streaming (for large files) |
| scala-json | JSON serialization/deserialization (for PDS upload) |
| org.scalafx / scalafx | GUI toolkit for the user interface (if using ScalaFX) |

### **1.2 Resource Management**

The application will use the src/main/resources directory to store internal data.

* **Reference Placeholder:** The core application logic will reference a designated path for the T2T-CHM13v2.0 reference files. The download logic must ensure these files (.fa, .fai, .dict) are placed in the application's local resource folder upon first use.

## **🔬 2\. High-Level Analysis Algorithms**

The core analysis is performed by the CallableLoci module, which replaces the former Ruby processing pipeline.

### **2.1 Algorithm: Callable Loci Analysis (GATK Integration)**

The Scala application performs the following high-level steps:

1. **Input Validation:** Verify the input BAM/CRAM file exists and is accessible.  
2. **Reference Check:** Check the local application resource folder for the required T2T-CHM13v2.0 reference files (.fa, .fai, .dict).  
   * **If missing:** Initiate a controlled, buffered download of the reference files from a known repository (e.g., NCBI, Broad).  
   * **If indices missing:** Use the HTSJDK API (CreateSequenceDictionary, Faidx generation) to generate the required indices locally.  
3. **GATK Execution (CallableLoci):** Execute the GATK CallableLoci tool programmatically, passing the input BAM and the local reference path.  
   * **Per-Contig Execution:** Run the analysis sequentially for each contig to manage memory and produce intermediate outputs.  
4. **Binning & Aggregation:** Read the output BED files from GATK.  
   * Iterate over the BED intervals and map them to fixed-size genomic bins (e.g., 10kb or 10Mb).  
   * For each bin, aggregate the status (Callable, Poor Coverage, No Coverage) and calculate the total base count for the histogram.  
5. **Data Structuring:** Compile the final data into two outputs:  
   * **a) SVG Generation:** Generate the per-contig SVG files using the dynamic width and fixed height logic derived from the prototype.  
   * **b) JSON Summary:** Generate a final JSON file containing the aggregated summary statistics and required metadata (pds\_user\_id, platform\_source).

### **2.2 Algorithm: JSON Summary for PDS Upload**

This is a critical algorithm for the Edge computing service:

1. **Metadata Collection:**  
   * platform\_source: Read from the BAM header or prompt the user.  
   * pds\_user\_id: Retrieved from a local configuration or a secure API call to decoding-us.com.  
2. **Summary Calculation:** Calculate genome-wide totals:  
   * Total Bases in Genome (from .dict file).  
   * Total Callable Bases (sum of all contigs).  
   * Total Poor Coverage Bases (sum of all contigs).  
   * Overall Callable Percentage (Total Callable / Total Bases).  
3. **JSON Serialization:** Serialize the aggregated data structure into the final JSON payload.  
4. **Secure Transmission:** Use a dedicated API client to encrypt and upload the JSON payload to the authenticated user's PDS data vault.

## **🖥️ 3\. UI/UX Screen Mockups (ScalaFX/Tauri Context)**

The application will use a simple, modern GUI to manage the entire process, abstracting the command line completely. The goal is a three-step process: Input, Analysis, and Report.

### **Screen A: Input & Welcome**

This screen focuses on simplicity and clear action for the user.

### **Screen B: Analysis Progress**

This screen provides continuous feedback, preventing user anxiety during long computation times.

### **Screen C: Results & PDS Opt-in**

The final screen presents actionable data and the crucial privacy decision for PDS contribution. The SVG plots are displayed in an integrated viewer.
