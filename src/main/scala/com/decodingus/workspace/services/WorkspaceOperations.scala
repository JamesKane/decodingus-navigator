package com.decodingus.workspace.services

import com.decodingus.workspace.model.*
import com.decodingus.workspace.{WorkspaceService, WorkspaceState}

/**
 * Handles all CRUD operations for workspace entities (Biosamples, Projects, SequenceRuns, Alignments).
 *
 * This service is stateless - it operates on the WorkspaceState provided and returns updated state.
 * The caller (WorkbenchViewModel) is responsible for persisting changes.
 */
class WorkspaceOperations {

  // --- Subject (Biosample) Operations ---

  def addSubject(state: WorkspaceState, newBiosample: Biosample, userDid: Option[String]): (WorkspaceState, Biosample) = {
    val enrichedBiosample = userDid match {
      case Some(did) =>
        val atUri = s"at://$did/com.decodingus.atmosphere.biosample/${newBiosample.sampleAccession}"
        newBiosample.copy(atUri = Some(atUri))
      case None =>
        newBiosample
    }

    val updatedSamples = state.workspace.main.samples :+ enrichedBiosample
    val updatedWorkspace = state.workspace.copy(main = state.workspace.main.copy(samples = updatedSamples))
    (state.copy(workspace = updatedWorkspace), enrichedBiosample)
  }

  def updateSubject(state: WorkspaceState, updatedBiosample: Biosample): WorkspaceState = {
    val updatedSamples = state.workspace.main.samples.map { sample =>
      if (sample.sampleAccession == updatedBiosample.sampleAccession) {
        updatedBiosample.copy(meta = sample.meta.updated("edit"))
      } else sample
    }
    state.copy(workspace = state.workspace.copy(main = state.workspace.main.copy(samples = updatedSamples)))
  }

  def updateSubjectDirect(state: WorkspaceState, updatedBiosample: Biosample): WorkspaceState = {
    val updatedSamples = state.workspace.main.samples.map { sample =>
      if (sample.sampleAccession == updatedBiosample.sampleAccession) updatedBiosample
      else sample
    }
    state.copy(workspace = state.workspace.copy(main = state.workspace.main.copy(samples = updatedSamples)))
  }

  def deleteSubject(state: WorkspaceState, sampleAccession: String): WorkspaceState = {
    val updatedSamples = state.workspace.main.samples.filterNot(_.sampleAccession == sampleAccession)
    state.copy(workspace = state.workspace.copy(main = state.workspace.main.copy(samples = updatedSamples)))
  }

  def findSubject(state: WorkspaceState, sampleAccession: String): Option[Biosample] = {
    state.workspace.main.samples.find(_.sampleAccession == sampleAccession)
  }

  // --- Project Operations ---

  def addProject(state: WorkspaceState, newProject: Project, userDid: Option[String]): (WorkspaceState, Project) = {
    val enrichedProject = userDid match {
      case Some(did) =>
        val rkey = java.util.UUID.randomUUID().toString
        val atUri = s"at://$did/com.decodingus.atmosphere.project/$rkey"
        newProject.copy(atUri = Some(atUri))
      case None =>
        newProject
    }

    val updatedProjects = state.workspace.main.projects :+ enrichedProject
    val updatedWorkspace = state.workspace.copy(main = state.workspace.main.copy(projects = updatedProjects))
    (state.copy(workspace = updatedWorkspace), enrichedProject)
  }

  def updateProject(state: WorkspaceState, updatedProject: Project): WorkspaceState = {
    val updatedProjects = state.workspace.main.projects.map { project =>
      if (project.projectName == updatedProject.projectName) {
        updatedProject.copy(meta = project.meta.updated("edit"))
      } else project
    }
    state.copy(workspace = state.workspace.copy(main = state.workspace.main.copy(projects = updatedProjects)))
  }

  def updateProjectDirect(state: WorkspaceState, updatedProject: Project): WorkspaceState = {
    val updatedProjects = state.workspace.main.projects.map { project =>
      if (project.projectName == updatedProject.projectName) updatedProject
      else project
    }
    state.copy(workspace = state.workspace.copy(main = state.workspace.main.copy(projects = updatedProjects)))
  }

  def deleteProject(state: WorkspaceState, projectName: String): WorkspaceState = {
    val updatedProjects = state.workspace.main.projects.filterNot(_.projectName == projectName)
    state.copy(workspace = state.workspace.copy(main = state.workspace.main.copy(projects = updatedProjects)))
  }

  def findProject(state: WorkspaceState, projectName: String): Option[Project] = {
    state.workspace.main.projects.find(_.projectName == projectName)
  }

  def addSubjectToProject(state: WorkspaceState, projectName: String, sampleAccession: String): Either[String, WorkspaceState] = {
    findProject(state, projectName) match {
      case Some(project) =>
        if (project.memberRefs.contains(sampleAccession)) {
          Left(s"Subject $sampleAccession already in project $projectName")
        } else {
          val updatedProject = project.copy(
            memberRefs = project.memberRefs :+ sampleAccession,
            meta = project.meta.updated("memberRefs")
          )
          Right(updateProjectDirect(state, updatedProject))
        }
      case None =>
        Left(s"Project $projectName not found")
    }
  }

  def removeSubjectFromProject(state: WorkspaceState, projectName: String, sampleAccession: String): Either[String, WorkspaceState] = {
    findProject(state, projectName) match {
      case Some(project) =>
        if (!project.memberRefs.contains(sampleAccession)) {
          Left(s"Subject $sampleAccession not in project $projectName")
        } else {
          val updatedProject = project.copy(
            memberRefs = project.memberRefs.filterNot(_ == sampleAccession),
            meta = project.meta.updated("memberRefs")
          )
          Right(updateProjectDirect(state, updatedProject))
        }
      case None =>
        Left(s"Project $projectName not found")
    }
  }

  def getProjectMembers(state: WorkspaceState, projectName: String): List[Biosample] = {
    findProject(state, projectName) match {
      case Some(project) =>
        project.memberRefs.flatMap(accession => findSubject(state, accession))
      case None =>
        List.empty
    }
  }

  def getNonProjectMembers(state: WorkspaceState, projectName: String): List[Biosample] = {
    findProject(state, projectName) match {
      case Some(project) =>
        state.workspace.main.samples.filterNot(s => project.memberRefs.contains(s.sampleAccession))
      case None =>
        state.workspace.main.samples
    }
  }

  // --- SequenceRun Operations ---

  /**
   * Creates a new SequenceRun from a FileInfo.
   * Returns updated state, the new run, and its index, or an error message.
   */
  def addSequenceRunFromFile(
    state: WorkspaceState,
    sampleAccession: String,
    fileInfo: FileInfo
  ): Either[String, (WorkspaceState, SequenceRun, Int)] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = state.workspace.main.getSequenceRunsForBiosample(subject)
        val existingChecksums = sequenceRuns.flatMap(_.files.flatMap(_.checksum)).toSet

        if (fileInfo.checksum.exists(existingChecksums.contains)) {
          Left(s"Duplicate file detected: ${fileInfo.fileName}")
        } else {
          val seqRunUri = s"local:sequencerun:${subject.sampleAccession}:${java.util.UUID.randomUUID().toString.take(8)}"

          val newSequenceRun = SequenceRun(
            atUri = Some(seqRunUri),
            meta = RecordMeta.initial,
            biosampleRef = subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}"),
            platformName = "Unknown",
            instrumentModel = None,
            testType = "Unknown",
            libraryLayout = None,
            totalReads = None,
            readLength = None,
            meanInsertSize = None,
            files = List(fileInfo),
            alignmentRefs = List.empty
          )

          val newIndex = sequenceRuns.size
          val updatedSequenceRuns = state.workspace.main.sequenceRuns :+ newSequenceRun
          val updatedSubject = subject.copy(
            sequenceRunRefs = subject.sequenceRunRefs :+ seqRunUri,
            meta = subject.meta.updated("sequenceRunRefs")
          )
          val updatedSamples = state.workspace.main.samples.map { s =>
            if (s.sampleAccession == sampleAccession) updatedSubject else s
          }
          val updatedContent = state.workspace.main.copy(
            samples = updatedSamples,
            sequenceRuns = updatedSequenceRuns
          )
          val newState = state.copy(workspace = state.workspace.copy(main = updatedContent))
          Right((newState, newSequenceRun, newIndex))
        }

      case None =>
        Left(s"Subject $sampleAccession not found")
    }
  }

  def getSequenceRun(state: WorkspaceState, sampleAccession: String, index: Int): Option[SequenceRun] = {
    findSubject(state, sampleAccession).flatMap { subject =>
      state.workspace.main.getSequenceRunsForBiosample(subject).lift(index)
    }
  }

  def getSequenceRunsForSubject(state: WorkspaceState, sampleAccession: String): List[SequenceRun] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) => state.workspace.main.getSequenceRunsForBiosample(subject)
      case None => List.empty
    }
  }

  def removeSequenceRun(state: WorkspaceState, sampleAccession: String, index: Int): Either[String, WorkspaceState] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = state.workspace.main.getSequenceRunsForBiosample(subject)
        if (index >= 0 && index < sequenceRuns.size) {
          val seqRunToRemove = sequenceRuns(index)

          val updatedSequenceRuns = state.workspace.main.sequenceRuns.filterNot(_.atUri == seqRunToRemove.atUri)
          val updatedAlignments = state.workspace.main.alignments.filterNot { align =>
            seqRunToRemove.atUri.exists(uri => align.sequenceRunRef == uri)
          }
          val updatedSubject = subject.copy(
            sequenceRunRefs = subject.sequenceRunRefs.filterNot(ref => seqRunToRemove.atUri.contains(ref)),
            meta = subject.meta.updated("sequenceRunRefs")
          )
          val updatedSamples = state.workspace.main.samples.map { s =>
            if (s.sampleAccession == sampleAccession) updatedSubject else s
          }
          val updatedContent = state.workspace.main.copy(
            samples = updatedSamples,
            sequenceRuns = updatedSequenceRuns,
            alignments = updatedAlignments
          )
          Right(state.copy(workspace = state.workspace.copy(main = updatedContent)))
        } else {
          Left(s"Index $index out of bounds")
        }
      case None =>
        Left(s"Subject $sampleAccession not found")
    }
  }

  def updateSequenceRun(state: WorkspaceState, sampleAccession: String, index: Int, updatedRun: SequenceRun): Either[String, WorkspaceState] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = state.workspace.main.getSequenceRunsForBiosample(subject)
        if (index >= 0 && index < sequenceRuns.size) {
          val originalRun = sequenceRuns(index)
          val runWithUpdatedMeta = updatedRun.copy(
            atUri = originalRun.atUri,
            meta = originalRun.meta.updated("edit")
          )
          val updatedSequenceRuns = state.workspace.main.sequenceRuns.map { sr =>
            if (sr.atUri == originalRun.atUri) runWithUpdatedMeta else sr
          }
          val updatedContent = state.workspace.main.copy(sequenceRuns = updatedSequenceRuns)
          Right(state.copy(workspace = state.workspace.copy(main = updatedContent)))
        } else {
          Left(s"Index $index out of bounds")
        }
      case None =>
        Left(s"Subject $sampleAccession not found")
    }
  }

  /**
   * Updates a sequence run by its URI directly.
   */
  def updateSequenceRunByUri(state: WorkspaceState, updatedRun: SequenceRun): WorkspaceState = {
    val updatedSequenceRuns = state.workspace.main.sequenceRuns.map { sr =>
      if (sr.atUri == updatedRun.atUri) updatedRun else sr
    }
    val updatedContent = state.workspace.main.copy(sequenceRuns = updatedSequenceRuns)
    state.copy(workspace = state.workspace.copy(main = updatedContent))
  }

  // --- Alignment Operations ---

  def addAlignment(state: WorkspaceState, alignment: Alignment): WorkspaceState = {
    val updatedAlignments = state.workspace.main.alignments :+ alignment
    val updatedContent = state.workspace.main.copy(alignments = updatedAlignments)
    state.copy(workspace = state.workspace.copy(main = updatedContent))
  }

  def updateAlignment(state: WorkspaceState, updatedAlignment: Alignment): WorkspaceState = {
    val updatedAlignments = state.workspace.main.alignments.map { a =>
      if (a.atUri == updatedAlignment.atUri) updatedAlignment else a
    }
    val updatedContent = state.workspace.main.copy(alignments = updatedAlignments)
    state.copy(workspace = state.workspace.copy(main = updatedContent))
  }

  def findAlignmentByUri(state: WorkspaceState, uri: String): Option[Alignment] = {
    state.workspace.main.alignments.find(_.atUri.contains(uri))
  }

  def getAlignmentsForSequenceRun(state: WorkspaceState, seqRun: SequenceRun): List[Alignment] = {
    state.workspace.main.getAlignmentsForSequenceRun(seqRun)
  }

  // --- AT URI Backfill ---

  def backfillAtUris(state: WorkspaceState, did: String): WorkspaceState = {
    var updated = false

    val updatedSamples = state.workspace.main.samples.map { sample =>
      if (sample.atUri.isEmpty) {
        updated = true
        val atUri = s"at://$did/com.decodingus.atmosphere.biosample/${sample.sampleAccession}"
        sample.copy(atUri = Some(atUri))
      } else sample
    }

    val updatedProjects = state.workspace.main.projects.map { project =>
      if (project.atUri.isEmpty) {
        updated = true
        val rkey = java.util.UUID.randomUUID().toString
        val atUri = s"at://$did/com.decodingus.atmosphere.project/$rkey"
        project.copy(atUri = Some(atUri))
      } else project
    }

    if (updated) {
      state.copy(workspace = state.workspace.copy(
        main = state.workspace.main.copy(samples = updatedSamples, projects = updatedProjects)
      ))
    } else {
      state
    }
  }

  // --- STR Profile Operations ---

  def addStrProfile(state: WorkspaceState, sampleAccession: String, profile: StrProfile): Either[String, (WorkspaceState, String)] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val strProfileUri = s"local:strprofile:${subject.sampleAccession}:${java.util.UUID.randomUUID().toString.take(8)}"
        val enrichedProfile = profile.copy(
          atUri = Some(strProfileUri),
          biosampleRef = subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}")
        )

        val updatedStrProfiles = state.workspace.main.strProfiles :+ enrichedProfile
        val updatedSubject = subject.copy(
          strProfileRefs = subject.strProfileRefs :+ strProfileUri,
          meta = subject.meta.updated("strProfileRefs")
        )
        val updatedSamples = state.workspace.main.samples.map { s =>
          if (s.sampleAccession == sampleAccession) updatedSubject else s
        }
        val updatedContent = state.workspace.main.copy(
          samples = updatedSamples,
          strProfiles = updatedStrProfiles
        )
        Right((state.copy(workspace = state.workspace.copy(main = updatedContent)), strProfileUri))

      case None =>
        Left(s"Subject not found: $sampleAccession")
    }
  }

  def getStrProfilesForBiosample(state: WorkspaceState, sampleAccession: String): List[StrProfile] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val byRefs = subject.strProfileRefs.flatMap { ref =>
          state.workspace.main.strProfiles.find(_.atUri.contains(ref))
        }
        if (byRefs.nonEmpty) byRefs
        else {
          val biosampleUri = subject.atUri.getOrElse(s"local:biosample:$sampleAccession")
          state.workspace.main.strProfiles.filter(_.biosampleRef == biosampleUri)
        }
      case None =>
        List.empty
    }
  }

  def getAllStrProfiles(state: WorkspaceState): List[StrProfile] = {
    state.workspace.main.strProfiles
  }

  def updateStrProfile(state: WorkspaceState, profileUri: String, updatedProfile: StrProfile): Either[String, WorkspaceState] = {
    state.workspace.main.strProfiles.find(_.atUri.contains(profileUri)) match {
      case Some(existing) =>
        val withUpdatedMeta = updatedProfile.copy(
          atUri = existing.atUri,
          meta = existing.meta.updated("edit")
        )
        val updatedStrProfiles = state.workspace.main.strProfiles.map { p =>
          if (p.atUri.contains(profileUri)) withUpdatedMeta else p
        }
        val updatedContent = state.workspace.main.copy(strProfiles = updatedStrProfiles)
        Right(state.copy(workspace = state.workspace.copy(main = updatedContent)))

      case None =>
        Left(s"STR profile not found: $profileUri")
    }
  }

  def deleteStrProfile(state: WorkspaceState, sampleAccession: String, profileUri: String): Either[String, WorkspaceState] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val updatedStrProfiles = state.workspace.main.strProfiles.filterNot(_.atUri.contains(profileUri))
        val updatedSubject = subject.copy(
          strProfileRefs = subject.strProfileRefs.filterNot(_ == profileUri),
          meta = subject.meta.updated("strProfileRefs")
        )
        val updatedSamples = state.workspace.main.samples.map { s =>
          if (s.sampleAccession == sampleAccession) updatedSubject else s
        }
        val updatedContent = state.workspace.main.copy(
          samples = updatedSamples,
          strProfiles = updatedStrProfiles
        )
        Right(state.copy(workspace = state.workspace.copy(main = updatedContent)))

      case None =>
        Left(s"Subject not found: $sampleAccession")
    }
  }

  // --- Chip Profile Operations ---

  def addChipProfile(state: WorkspaceState, sampleAccession: String, profile: ChipProfile): Either[String, (WorkspaceState, String)] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val chipProfileUri = profile.atUri.getOrElse(
          s"local:chipprofile:${subject.sampleAccession}:${java.util.UUID.randomUUID().toString.take(8)}"
        )
        val enrichedProfile = profile.copy(
          atUri = Some(chipProfileUri),
          biosampleRef = subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}")
        )

        val updatedChipProfiles = state.workspace.main.chipProfiles :+ enrichedProfile
        val updatedSubject = subject.copy(
          genotypeRefs = subject.genotypeRefs :+ chipProfileUri,
          meta = subject.meta.updated("genotypeRefs")
        )
        val updatedSamples = state.workspace.main.samples.map { s =>
          if (s.sampleAccession == sampleAccession) updatedSubject else s
        }
        val updatedContent = state.workspace.main.copy(
          samples = updatedSamples,
          chipProfiles = updatedChipProfiles
        )
        Right((state.copy(workspace = state.workspace.copy(main = updatedContent)), chipProfileUri))

      case None =>
        Left(s"Subject not found: $sampleAccession")
    }
  }

  def getChipProfilesForBiosample(state: WorkspaceState, sampleAccession: String): List[ChipProfile] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val byRefs = subject.genotypeRefs.flatMap { ref =>
          state.workspace.main.chipProfiles.find(_.atUri.contains(ref))
        }
        if (byRefs.nonEmpty) byRefs
        else {
          val biosampleUri = subject.atUri.getOrElse(s"local:biosample:$sampleAccession")
          state.workspace.main.chipProfiles.filter(_.biosampleRef == biosampleUri)
        }
      case None =>
        List.empty
    }
  }

  def deleteChipProfile(state: WorkspaceState, sampleAccession: String, profileUri: String): Either[String, WorkspaceState] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val updatedChipProfiles = state.workspace.main.chipProfiles.filterNot(_.atUri.contains(profileUri))
        val updatedSubject = subject.copy(
          genotypeRefs = subject.genotypeRefs.filterNot(_ == profileUri),
          meta = subject.meta.updated("genotypeRefs")
        )
        val updatedSamples = state.workspace.main.samples.map { s =>
          if (s.sampleAccession == sampleAccession) updatedSubject else s
        }
        val updatedContent = state.workspace.main.copy(
          samples = updatedSamples,
          chipProfiles = updatedChipProfiles
        )
        Right(state.copy(workspace = state.workspace.copy(main = updatedContent)))

      case None =>
        Left(s"Subject not found: $sampleAccession")
    }
  }

  // --- Utility ---

  def getExistingChecksums(state: WorkspaceState, sampleAccession: String): Set[String] = {
    findSubject(state, sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = state.workspace.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.flatMap(_.files.flatMap(_.checksum)).toSet
      case None =>
        Set.empty
    }
  }
}
