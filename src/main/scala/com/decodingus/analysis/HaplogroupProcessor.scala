package com.decodingus.analysis

import com.decodingus.config.FeatureToggles
import com.decodingus.haplogroup.caller.{GatkHaplotypeCallerProcessor, TwoPassCallerResult}
import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupResult, Locus}
import com.decodingus.haplogroup.report.HaplogroupReportWriter
import com.decodingus.haplogroup.scoring.HaplogroupScorer
import com.decodingus.haplogroup.tree.{TreeCache, TreeProvider, TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider, NamedVariantCache}
import com.decodingus.liftover.LiftoverProcessor
import com.decodingus.model.LibraryStats
import com.decodingus.refgenome.{LiftoverGateway, MultiContigReferenceQuerier, ReferenceGateway, ReferenceQuerier, StrAnnotator}
import htsjdk.variant.vcf.VCFFileReader

import java.io.{File, PrintWriter}
import java.nio.file.{Files, Path}
import scala.jdk.CollectionConverters.*
import scala.util.Using

/**
 * A private/novel variant not found in the haplogroup tree.
 */
case class PrivateVariant(
  contig: String,
  position: Long,
  ref: String,
  alt: String,
  quality: Option[Double]
)

class HaplogroupProcessor {

  private val standardContigOrder: Map[String, Int] = (1 to 22).map(i => s"chr$i" -> i).toMap ++
    Map("chrX" -> 23, "chrY" -> 24, "chrM" -> 25)

  private val ARTIFACT_SUBDIR_NAME = "haplogroup"

  /**
   * Analyze a BAM/CRAM file for haplogroup assignment.
   *
   * @param bamPath Path to the BAM/CRAM file
   * @param libraryStats Library statistics from initial analysis
   * @param treeType Y-DNA or MT-DNA tree type
   * @param treeProviderType Tree data provider (FTDNA or DecodingUs)
   * @param onProgress Progress callback
   * @param artifactContext Optional context for organizing output artifacts by subject/run/alignment
   */
  def analyze(
               bamPath: String,
               libraryStats: LibraryStats,
               treeType: TreeType,
               treeProviderType: TreeProviderType,
               onProgress: (String, Double, Double) => Unit,
               artifactContext: Option[ArtifactContext] = None
             ): Either[String, List[HaplogroupResult]] = {

    onProgress("Loading haplogroup tree...", 0.0, 1.0)
    val treeProvider: TreeProvider = treeProviderType match {
      case TreeProviderType.FTDNA => new FtdnaTreeProvider(treeType)
      case TreeProviderType.DECODINGUS => new DecodingUsTreeProvider(treeType)
    }
    val treeCache = new TreeCache()

    treeProvider.loadTree(libraryStats.referenceBuild).flatMap { tree =>
      val allLoci = collectAllLoci(tree).distinct
      // The 'contig' variable is no longer directly used here for createVcfAllelesFile,
      // as Locus objects now carry their own contig information.
      // However, it's still used for liftover and other contig-specific operations.
      val primaryContig = if (treeType == TreeType.YDNA) "chrY" else "chrM"
      val outputPrefix = if (treeType == TreeType.YDNA) "ydna" else "mtdna"

      val referenceBuild = libraryStats.referenceBuild
      val treeSourceBuild = if (treeProvider.supportedBuilds.contains(referenceBuild)) {
        referenceBuild
      } else {
        treeProvider.sourceBuild
      }

      // Check if we can optimize by using reference genome's known haplogroup
      // Only applies to Y-DNA currently (MT-DNA reference sources are less characterized)
      val (pass1Loci, referenceHaplogroup) = if (treeType == TreeType.YDNA) {
        FeatureToggles.referenceHaplogroups.getHaplogroups(referenceBuild) match {
          case Some(haplogroupNames) =>
            // Try to find the reference haplogroup in the tree (using any of the name variants)
            val foundHaplogroup = haplogroupNames.flatMap(name => findHaplogroupByName(tree, name)).headOption
            foundHaplogroup match {
              case Some(refHg) =>
                // Collect only loci on the path from root to reference haplogroup
                val pathLoci = collectPathLoci(tree, refHg.name)
                onProgress(s"Optimizing: Reference is ${refHg.name}, using ${pathLoci.size} path positions for pass 1...", 0.02, 1.0)
                (pathLoci, Some(refHg.name))
              case None =>
                onProgress(s"Reference haplogroup not found in tree, using full tree (${allLoci.size} positions)...", 0.02, 1.0)
                (allLoci, None)
            }
          case None =>
            onProgress(s"No known haplogroup for $referenceBuild, using full tree (${allLoci.size} positions)...", 0.02, 1.0)
            (allLoci, None)
        }
      } else {
        // MT-DNA: use all loci for now
        (allLoci, None)
      }

      val referenceGateway = new ReferenceGateway((_, _) => {})

      referenceGateway.resolve(treeSourceBuild).flatMap { treeRefPath =>
        // Cache key includes reference haplogroup for path-optimized VCFs
        val cacheKeySuffix = referenceHaplogroup.map(h => s"-path-$h").getOrElse("")
        val primaryCacheKey = s"${treeProvider.cachePrefix}$cacheKeySuffix"
        val cachedSitesVcf = treeCache.getSitesVcfPath(primaryCacheKey, treeSourceBuild)

        // Check for existing VCFs - try current provider's cache first, then check alternatives
        val initialAllelesVcf = if (treeCache.isSitesVcfValid(primaryCacheKey, treeSourceBuild)) {
          onProgress("Using cached sites VCF...", 0.05, 1.0)
          cachedSitesVcf
        } else {
          // Check for alternative cached VCFs (e.g., from other tree providers)
          // This allows reusing VCFs if someone switches providers
          val alternativePrefixes = treeType match {
            case TreeType.YDNA => List("ftdna-ytree", "decodingus-ytree")
            case TreeType.MTDNA => List("ftdna-mttree", "decodingus-mttree")
          }
          val alternativeVcf = alternativePrefixes
            .filterNot(_ == treeProvider.cachePrefix) // Skip current provider
            .map(prefix => s"$prefix$cacheKeySuffix")
            .map(key => treeCache.getSitesVcfPath(key, treeSourceBuild))
            .find(_.exists())

          alternativeVcf match {
            case Some(existingVcf) =>
              onProgress(s"Found existing sites VCF from alternative provider...", 0.05, 1.0)
              existingVcf
            case None =>
              onProgress(s"Creating sites VCF (${pass1Loci.size} positions)...", 0.05, 1.0)
              createVcfAllelesFile(pass1Loci, treeRefPath.toString, treeType, Some(cachedSitesVcf))
          }
        }

        // Define artifact directory early so it can be used for saving intermediate files
        val artifactDir = artifactContext.map(_.getSubdir(ARTIFACT_SUBDIR_NAME))

        val (allelesForCalling, performReverseLiftover) = if (referenceBuild == treeSourceBuild) {
          onProgress("Reference builds match.", 0.1, 1.0)
          (Right(initialAllelesVcf), false)
        } else {
          // Check for cached lifted VCF first (shared across all samples with same tree/build combo)
          val cachedLiftedVcf = treeCache.getLiftedSitesVcfPath(primaryCacheKey, treeSourceBuild, referenceBuild)
          if (treeCache.isLiftedSitesVcfValid(primaryCacheKey, treeSourceBuild, referenceBuild)) {
            onProgress(s"Using cached lifted sites VCF for $referenceBuild...", 0.1, 1.0)
            println(s"[HaplogroupProcessor] Using cached lifted sites VCF: $cachedLiftedVcf")
            (Right(cachedLiftedVcf), true)
          } else {
            onProgress(s"Reference mismatch: tree is $treeSourceBuild, BAM/CRAM is $referenceBuild. Performing liftover...", 0.1, 1.0)
            // Note: The contig parameter here for liftover still refers to the primary contig for the tree type.
            // Filter output to only keep expected contig (chrY or chrM) to exclude NUMT mappings
            val lifted = performLiftover(initialAllelesVcf, primaryContig, treeSourceBuild, referenceBuild, onProgress, filterOutput = true)
            // Cache the lifted VCF for reuse by other samples
            lifted.foreach { liftedVcf =>
              java.nio.file.Files.copy(liftedVcf.toPath, cachedLiftedVcf.toPath, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
              println(s"[HaplogroupProcessor] Cached lifted sites VCF to $cachedLiftedVcf")
            }
            // Also save a copy to artifact directory for debugging/auditing
            lifted.foreach { liftedVcf =>
              artifactDir.foreach { dir =>
                val savedPath = dir.resolve(s"${outputPrefix}_lifted_alleles_${referenceBuild}.vcf")
                java.nio.file.Files.copy(liftedVcf.toPath, savedPath, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
                println(s"[HaplogroupProcessor] Saved lifted alleles VCF to $savedPath")
              }
            }
            (lifted, true)
          }
        }

        allelesForCalling.flatMap { vcf =>
          referenceGateway.resolve(referenceBuild).flatMap { referencePath =>
            val caller = new GatkHaplotypeCallerProcessor()

            // Two-pass calling: tree sites first, then private variants
            caller.callTwoPass(
              bamPath,
              referencePath.toString,
              vcf,
              (msg, done, total) => onProgress(msg, 0.2 + (done * 0.5), 1.0),
              artifactDir,
              Some(outputPrefix),
              Some(primaryContig) // Explicitly pass the primary contig (chrY or chrM)
            ).flatMap { twoPassResult =>
              val postGatkStart = System.currentTimeMillis()
              // Handle reverse liftover for tree sites VCF if needed
              // Filter to primaryContig to remove variants that mapped to unexpected contigs (e.g., chrY -> chrX in PAR)
              val finalTreeVcf = if (performReverseLiftover) {
                onProgress("Performing reverse liftover on tree sites...", 0.72, 1.0)
                val reverseLifted = performLiftover(twoPassResult.treeSitesVcf, primaryContig, referenceBuild, treeSourceBuild, onProgress, filterOutput = true)
                // Save reverse-lifted VCF (used for scoring) to artifact directory
                reverseLifted.foreach { liftedVcf =>
                  artifactDir.foreach { dir =>
                    val savedPath = dir.resolve(s"${outputPrefix}_calls_lifted_${treeSourceBuild}.vcf")
                    java.nio.file.Files.copy(liftedVcf.toPath, savedPath, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
                    println(s"[HaplogroupProcessor] Saved reverse-lifted calls VCF to $savedPath")
                  }
                }
                reverseLifted
              } else {
                Right(twoPassResult.treeSitesVcf)
              }

              // Also lift private variants VCF back to tree coordinates if needed
              val finalPrivateVariantsVcf = if (performReverseLiftover) {
                onProgress("Lifting private variants back to tree coordinates...", 0.76, 1.0)
                val liftedPrivate = performLiftover(twoPassResult.privateVariantsVcf, primaryContig, referenceBuild, treeSourceBuild, onProgress, filterOutput = true)
                // Save lifted private variants VCF
                liftedPrivate.foreach { liftedVcf =>
                  artifactDir.foreach { dir =>
                    val savedPath = dir.resolve(s"${outputPrefix}_private_variants_lifted_${treeSourceBuild}.vcf")
                    java.nio.file.Files.copy(liftedVcf.toPath, savedPath, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
                    println(s"[HaplogroupProcessor] Saved lifted private variants VCF to $savedPath")
                  }
                }
                liftedPrivate
              } else {
                Right(twoPassResult.privateVariantsVcf)
              }

              finalTreeVcf.flatMap { scoredVcf =>
                finalPrivateVariantsVcf.flatMap { privateVcf =>
                  onProgress("Scoring haplogroups...", 0.8, 1.0)
                  val scoringStart = System.currentTimeMillis()
                  // Merge calls from both VCFs:
                  // 1. Tree sites VCF (pass 1) - forced calls at reference-path positions
                  // 2. Private variants VCF (pass 2) - variant calls across full chromosome
                  // Pass 2 may contain calls at tree positions not in pass 1 (e.g., other branches)
                  // Note: Both VCFs are now in tree source coordinates for proper position matching
                  val allTreePositions = allLoci.map(_.position).toSet
                  val treeSiteCalls = parseVcf(scoredVcf)
                  println(s"[HaplogroupProcessor] Parsed tree sites VCF: ${treeSiteCalls.size} calls in ${System.currentTimeMillis() - scoringStart}ms")
                  val additionalTreeCalls = parseVcfAtPositions(privateVcf, allTreePositions)
                  println(s"[HaplogroupProcessor] Parsed additional calls in ${System.currentTimeMillis() - scoringStart}ms")
                  // Merge: pass 1 calls take precedence (they're force-called at exact positions)
                  val snpCalls = additionalTreeCalls ++ treeSiteCalls

                  val scorer = new HaplogroupScorer()
                  val scoreStart = System.currentTimeMillis()
                  val results = scorer.score(tree, snpCalls)
                  println(s"[HaplogroupProcessor] Scored ${results.size} haplogroups in ${System.currentTimeMillis() - scoreStart}ms")

                  // Identify private variants - only exclude positions on path to terminal haplogroup
                  // Positions on other branches could be legitimate private variants for undiscovered sub-clades
                  // Note: Use lifted VCF for position matching, but original VCF coordinates are kept for reporting
                  onProgress("Identifying private variants...", 0.85, 1.0)
                  val terminalHaplogroup = results.headOption.map(_.name).getOrElse("")
                  val pathPositions = collectPathPositions(tree, terminalHaplogroup)
                  // Use the original (non-lifted) private variants for the report - positions in BAM's reference
                  val privateVariants = parsePrivateVariants(twoPassResult.privateVariantsVcf, pathPositions)
                  println(s"[HaplogroupProcessor] Post-GATK processing completed in ${System.currentTimeMillis() - postGatkStart}ms")

                  // Load STR annotator for indel annotation (optional - don't fail if unavailable)
                  val strAnnotator = StrAnnotator.forBuild(referenceBuild) match {
                    case Right(annotator) =>
                      println(s"[HaplogroupProcessor] Loaded STR reference with ${annotator.regionCount} regions")
                      Some(annotator)
                    case Left(error) =>
                      println(s"[HaplogroupProcessor] STR annotation unavailable: $error")
                      None
                  }

                  // Write report to artifact directory if available
                  artifactDir.foreach { dir =>
                    onProgress("Writing haplogroup report...", 0.9, 1.0)

                    // Use named variant cache for Decoding Us provider to enrich reports with aliases
                    val variantCache = treeProviderType match {
                      case TreeProviderType.DECODINGUS =>
                        val cache = NamedVariantCache()
                        // Try to load silently - don't fail the report if cache unavailable
                        cache.ensureLoaded(msg => println(s"[HaplogroupProcessor] $msg")) match {
                          case Right(_) => Some(cache)
                          case Left(err) =>
                            println(s"[HaplogroupProcessor] Named variant cache unavailable: $err")
                            None
                        }
                      case _ => None
                    }

                    HaplogroupReportWriter.writeReport(
                      outputDir = dir.toFile,
                      treeType = treeType,
                      results = results,
                      tree = tree,
                      snpCalls = snpCalls,
                      sampleName = None,
                      privateVariants = Some(privateVariants),
                      treeProvider = Some(treeProviderType),
                      strAnnotator = strAnnotator,
                      sampleBuild = Some(referenceBuild),
                      treeBuild = Some(treeSourceBuild),
                      namedVariantCache = variantCache
                    )
                  }

                  onProgress("Analysis complete.", 1.0, 1.0)
                  Right(results)
                }
              }
            }
          }
        }
      }
    }
  }

  /**
   * Analyze using a cached whole-genome VCF instead of calling variants on the fly.
   * Uses GapAwareHaplogroupResolver to query the VCF and infer reference calls from callable loci.
   *
   * @param sampleAccession Sample accession for artifact lookup
   * @param runId Sequence run ID
   * @param alignmentId Alignment ID
   * @param referenceBuild Reference build of the alignment
   * @param treeType Y-DNA or MT-DNA tree type
   * @param treeProviderType Tree data provider
   * @param onProgress Progress callback
   */
  def analyzeFromCachedVcf(
    sampleAccession: String,
    runId: String,
    alignmentId: String,
    referenceBuild: String,
    treeType: TreeType,
    treeProviderType: TreeProviderType,
    onProgress: (String, Double, Double) => Unit
  ): Either[String, List[HaplogroupResult]] = {

    onProgress("Loading haplogroup tree...", 0.0, 1.0)
    val treeProvider: TreeProvider = treeProviderType match {
      case TreeProviderType.FTDNA => new FtdnaTreeProvider(treeType)
      case TreeProviderType.DECODINGUS => new DecodingUsTreeProvider(treeType)
    }

    treeProvider.loadTree(referenceBuild).flatMap { tree =>
      val allLoci = collectAllLoci(tree).distinct
      val primaryContig = if (treeType == TreeType.YDNA) "chrY" else "chrM"
      val outputPrefix = if (treeType == TreeType.YDNA) "ydna" else "mtdna"

      val treeSourceBuild = if (treeProvider.supportedBuilds.contains(referenceBuild)) {
        referenceBuild
      } else {
        treeProvider.sourceBuild
      }

      onProgress("Loading cached VCF and callable loci...", 0.1, 1.0)

      // Create gap-aware resolver
      GapAwareHaplogroupResolver.fromCache(sampleAccession, runId, alignmentId, referenceBuild) match {
        case Left(error) =>
          Left(s"Failed to load cached VCF: $error")

        case Right(resolver) =>
          onProgress("Querying tree positions from VCF...", 0.2, 1.0)

          // Query all tree positions
          val positions = allLoci.map { locus =>
            val contig = if (locus.contig.isEmpty) primaryContig else locus.contig
            (contig, locus.position, locus.ref)
          }
          val resolvedCalls = resolver.resolvePositions(positions)

          // Get resolution statistics
          val stats = resolver.getResolutionStats(resolvedCalls)
          onProgress(s"Resolved ${stats.fromVcf} from VCF, ${stats.inferredReference} inferred, ${stats.noCalls} no-calls", 0.4, 1.0)

          // Convert to snpCalls format expected by scorer (position -> allele)
          val snpCalls: Map[Long, String] = resolvedCalls.flatMap { case ((contig, pos), call) =>
            if (call.hasCall && contig == primaryContig) {
              Some(pos -> call.allele)
            } else {
              None
            }
          }

          onProgress("Scoring haplogroups...", 0.5, 1.0)

          // Score haplogroups
          val scorer = new HaplogroupScorer()
          val results = scorer.score(tree, snpCalls)

          onProgress("Generating report...", 0.8, 1.0)

          // Get artifact directory for report
          val outputDir = SubjectArtifactCache.getArtifactSubdir(sampleAccession, runId, alignmentId, ARTIFACT_SUBDIR_NAME)

          // Use named variant cache for Decoding Us provider
          val variantCache = treeProviderType match {
            case TreeProviderType.DECODINGUS =>
              val cache = NamedVariantCache()
              cache.ensureLoaded(_ => ()) match {
                case Right(_) => Some(cache)
                case Left(_) => None
              }
            case _ => None
          }

          HaplogroupReportWriter.writeReport(
            outputDir = outputDir.toFile,
            treeType = treeType,
            results = results,
            tree = tree,
            snpCalls = snpCalls,
            sampleName = None,
            privateVariants = None,  // Private variants not extracted in VCF-based flow
            treeProvider = Some(treeProviderType),
            strAnnotator = None,
            sampleBuild = Some(referenceBuild),
            treeBuild = Some(treeSourceBuild),
            namedVariantCache = variantCache
          )

          onProgress("Analysis complete.", 1.0, 1.0)
          Right(results)
      }
    }
  }

  /**
   * Parse private variants from VCF, excluding known tree positions.
   */
  private def parsePrivateVariants(vcfFile: File, treePositions: Set[Long]): List[PrivateVariant] = {
    val reader = new VCFFileReader(vcfFile, false)
    val variants = reader.iterator().asScala.flatMap { vc =>
      val pos = vc.getStart.toLong
      if (!treePositions.contains(pos)) {
        val genotype = vc.getGenotypes.get(0)
        val allele = genotype.getAlleles.get(0).getBaseString
        val ref = vc.getReference.getBaseString
        val qual = if (vc.hasLog10PError) Some(vc.getPhredScaledQual) else None
        Some(PrivateVariant(
          contig = vc.getContig,
          position = pos,
          ref = ref,
          alt = allele,
          quality = qual
        ))
      } else {
        None
      }
    }.toList
    reader.close()
    variants
  }

  /**
   * Perform liftover of a VCF file between reference builds.
   *
   * @param vcfFile Input VCF file
   * @param expectedContig Expected contig for filtering (chrY or chrM) - used to filter out
   *                       variants that map to unexpected contigs after liftover
   * @param fromBuild Source reference build
   * @param toBuild Target reference build
   * @param onProgress Progress callback
   * @param filterOutput If true, filter output to only include variants on expectedContig.
   *                     Should be true for reverse liftover (back to tree coordinates).
   */
  private def performLiftover(
                               vcfFile: File,
                               expectedContig: String,
                               fromBuild: String,
                               toBuild: String,
                               onProgress: (String, Double, Double) => Unit,
                               filterOutput: Boolean = false
                             ): Either[String, File] = {
    val liftoverGateway = new LiftoverGateway((_, _) => {})
    val referenceGateway = new ReferenceGateway((_, _) => {})

    val filterToContig = if (filterOutput) Some(expectedContig) else None

    for {
      chainFile <- liftoverGateway.resolve(fromBuild, toBuild)
      targetRef <- referenceGateway.resolve(toBuild)
      liftedVcf <- new LiftoverProcessor().liftoverVcf(vcfFile, chainFile, targetRef, (msg, done, total) => onProgress(msg, 0.2 + (done * 0.2), 1.0), filterToContig)
    } yield liftedVcf
  }

  private def createVcfAllelesFile(
    loci: List[Locus],
    referencePath: String,
    treeType: TreeType,
    outputFile: Option[File]
  ): File = {
    // Use provided output file or create temp file
    val vcfFile = outputFile match {
      case Some(file) =>
        // Ensure parent directory exists
        Option(file.getParentFile).foreach(_.mkdirs())
        file
      case None =>
        val tempFile = File.createTempFile("alleles", ".vcf")
        tempFile.deleteOnExit()
        tempFile
    }

    // Group loci by contig first, then by position within each contig
    val lociByContig = loci.groupBy(_.contig)
    val sortedContigs = lociByContig.keys.toList.sortBy(c => standardContigOrder.getOrElse(c, 999))

    Using.resource(new PrintWriter(vcfFile)) { writer =>
      writer.println("##fileformat=VCFv4.2")
      sortedContigs.foreach { c =>
        writer.println(s"##contig=<ID=$c>")
      }
      writer.println("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">")
      writer.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

      // Process one contig at a time - load reference once per contig
      sortedContigs.foreach { contig =>
        Using.resource(new ReferenceQuerier(referencePath, contig)) { refQuerier =>
          val contigLoci = lociByContig(contig)
          // Group by position to combine alternates at the same site
          val groupedByPosition = contigLoci.groupBy(_.position)
          // Filter out positions beyond contig bounds, then sort
          val sortedPositions = groupedByPosition.keys.toList
            .filter(refQuerier.isValidPosition)
            .sorted

          sortedPositions.foreach { position =>
            val lociAtPosition = groupedByPosition(position)
            // Safe to use get since we filtered valid positions above
            val refBase = refQuerier.getBase(position).get

            // Filter to valid SNPs only:
            // - Single base ref and alt (no indels)
            // - Only valid nucleotides A, C, G, T (no dashes, dots, or other characters)
            val validBases = Set("A", "C", "G", "T")
            val snpLoci = lociAtPosition.filter { l =>
              l.ref.length == 1 && l.alt.length == 1 &&
                validBases.contains(l.ref.toUpperCase) &&
                validBases.contains(l.alt.toUpperCase)
            }

            // Collect all valid alternates at this position
            val refAndAlts = snpLoci.flatMap { locus =>
              if (refBase.toString.equalsIgnoreCase(locus.ref)) {
                Some((locus.ref.toUpperCase, locus.alt.toUpperCase))
              } else if (refBase.toString.equalsIgnoreCase(locus.alt)) {
                Some((locus.alt.toUpperCase, locus.ref.toUpperCase))
              } else {
                // This locus is problematic, the reference doesn't match ANC or DER
                None
              }
            }

            if (refAndAlts.nonEmpty) {
              // All refs should be the same (the actual reference base)
              val ref = refAndAlts.head._1
              // Collect unique alternates
              val alts = refAndAlts.map(_._2).distinct.filterNot(_ == ref)

              if (alts.nonEmpty) {
                writer.println(s"$contig\t$position\t.\t$ref\t${alts.mkString(",")}\t.\t.\t.")
              }
            }
          }
        }
      }
    }
    vcfFile
  }

  /**
   * Recursively collect all loci from the haplogroup tree.
   */
  private def collectAllLoci(tree: List[Haplogroup]): List[Locus] = {
    tree.flatMap(collectLociFromHaplogroup)
  }

  private def collectLociFromHaplogroup(haplogroup: Haplogroup): List[Locus] = {
    haplogroup.loci ++ haplogroup.children.flatMap(collectLociFromHaplogroup)
  }

  /**
   * Find a haplogroup by name in the tree (case-insensitive).
   */
  private def findHaplogroupByName(tree: List[Haplogroup], name: String): Option[Haplogroup] = {
    def search(haplogroup: Haplogroup): Option[Haplogroup] = {
      if (haplogroup.name.equalsIgnoreCase(name)) {
        Some(haplogroup)
      } else {
        haplogroup.children.flatMap(search).headOption
      }
    }
    tree.flatMap(search).headOption
  }

  /**
   * Collect all loci along the path from root to the specified terminal haplogroup.
   * Used for optimized pass 1 calling when we know the reference genome's haplogroup.
   */
  private def collectPathLoci(tree: List[Haplogroup], terminalName: String): List[Locus] = {
    findPath(tree, terminalName).flatMap(_.loci).distinct
  }

  /**
   * Find the path from root to a haplogroup by name.
   */
  private def findPath(tree: List[Haplogroup], terminalName: String): List[Haplogroup] = {
    def findPathFromNode(haplogroup: Haplogroup): Option[List[Haplogroup]] = {
      if (haplogroup.name.equalsIgnoreCase(terminalName)) {
        Some(List(haplogroup))
      } else {
        haplogroup.children.flatMap(findPathFromNode).headOption.map(path => haplogroup :: path)
      }
    }
    tree.flatMap(findPathFromNode).headOption.getOrElse(List.empty)
  }

  /**
   * Collect all positions along the path from root to the specified terminal haplogroup.
   * Only these positions should be excluded from private variant detection.
   */
  private def collectPathPositions(tree: List[Haplogroup], terminalName: String): Set[Long] = {
    findPath(tree, terminalName).flatMap(_.loci.map(_.position)).toSet
  }

  private def parseVcf(vcfFile: File): Map[Long, String] = {
    val reader = new VCFFileReader(vcfFile, false)
    val snpCalls = reader.iterator().asScala.map {
      vc =>
        val pos = vc.getStart.toLong
        val genotype = vc.getGenotypes.get(0) // Assuming single sample VCF
        val allele = genotype.getAlleles.get(0).getBaseString
        pos -> allele
    }.toMap
    reader.close()
    snpCalls
  }

  /**
   * Parse VCF and return only calls at specified positions.
   * Used to extract tree-position calls from the private variants VCF.
   */
  private def parseVcfAtPositions(vcfFile: File, positions: Set[Long]): Map[Long, String] = {
    val reader = new VCFFileReader(vcfFile, false)
    val snpCalls = reader.iterator().asScala.flatMap { vc =>
      val pos = vc.getStart.toLong
      if (positions.contains(pos)) {
        val genotype = vc.getGenotypes.get(0)
        val allele = genotype.getAlleles.get(0).getBaseString
        Some(pos -> allele)
      } else {
        None
      }
    }.toMap
    reader.close()
    snpCalls
  }
}