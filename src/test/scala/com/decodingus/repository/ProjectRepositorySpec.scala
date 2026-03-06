package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import munit.FunSuite
import java.util.UUID

class ProjectRepositorySpec extends FunSuite with DatabaseTestSupport:

  val projectRepo = ProjectRepository()
  val biosampleRepo = BiosampleRepository()

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val entity = ProjectEntity.create(
        projectName = "Test Project",
        administratorDid = "did:plc:admin123",
        description = Some("A test project")
      )

      val saved = projectRepo.insert(entity)
      val found = projectRepo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.projectName, "Test Project")
      assertEquals(found.get.administratorDid, "did:plc:admin123")
      assertEquals(found.get.description, Some("A test project"))
      assertEquals(found.get.meta.syncStatus, SyncStatus.Local)
    }
  }

  testTransactor.test("findAll returns all projects") { case (db, tx) =>
    tx.readWrite {
      projectRepo.insert(ProjectEntity.create("Project A", "did:admin1"))
      projectRepo.insert(ProjectEntity.create("Project B", "did:admin1"))
      projectRepo.insert(ProjectEntity.create("Project C", "did:admin2"))

      val all = projectRepo.findAll()
      assertEquals(all.size, 3)
    }
  }

  testTransactor.test("update modifies entity correctly") { case (db, tx) =>
    tx.readWrite {
      val entity = projectRepo.insert(ProjectEntity.create("Original", "did:admin"))

      val updated = entity.copy(
        projectName = "Updated",
        description = Some("New description")
      )
      projectRepo.update(updated)

      val found = projectRepo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.projectName, "Updated")
      assertEquals(found.get.description, Some("New description"))
      assertEquals(found.get.meta.version, 2)
    }
  }

  testTransactor.test("delete removes project and memberships via cascade") { case (db, tx) =>
    tx.readWrite {
      // Create project with members
      val project = projectRepo.insert(ProjectEntity.create("DeleteMe", "did:admin"))
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      projectRepo.addMember(project.id, biosample.id)
      assert(projectRepo.isMember(project.id, biosample.id))

      // Delete project
      projectRepo.delete(project.id)

      // Project should be gone
      assertEquals(projectRepo.findById(project.id), None)

      // Membership should be cascade deleted
      assertEquals(projectRepo.getMemberIds(project.id), List.empty)
    }
  }

  testTransactor.test("findByName finds project by unique name") { case (db, tx) =>
    tx.readWrite {
      projectRepo.insert(ProjectEntity.create("UniqueProject", "did:admin"))

      val found = projectRepo.findByName("UniqueProject")
      assert(found.isDefined)
      assertEquals(found.get.administratorDid, "did:admin")

      assertEquals(projectRepo.findByName("NonExistent"), None)
    }
  }

  testTransactor.test("findByAdministrator returns projects for admin") { case (db, tx) =>
    tx.readWrite {
      projectRepo.insert(ProjectEntity.create("Admin1 Project A", "did:admin1"))
      projectRepo.insert(ProjectEntity.create("Admin1 Project B", "did:admin1"))
      projectRepo.insert(ProjectEntity.create("Admin2 Project", "did:admin2"))

      val admin1Projects = projectRepo.findByAdministrator("did:admin1")
      assertEquals(admin1Projects.size, 2)

      val admin2Projects = projectRepo.findByAdministrator("did:admin2")
      assertEquals(admin2Projects.size, 1)
    }
  }

  testTransactor.test("searchByName finds by prefix") { case (db, tx) =>
    tx.readWrite {
      projectRepo.insert(ProjectEntity.create("Alpha Project", "did:admin"))
      projectRepo.insert(ProjectEntity.create("Alpha Research", "did:admin"))
      projectRepo.insert(ProjectEntity.create("Beta Study", "did:admin"))

      val alphaResults = projectRepo.searchByName("Alpha")
      assertEquals(alphaResults.size, 2)

      val betaResults = projectRepo.searchByName("Beta")
      assertEquals(betaResults.size, 1)
    }
  }

  // ============================================
  // Membership Tests
  // ============================================

  testTransactor.test("addMember adds biosample to project") { case (db, tx) =>
    tx.readWrite {
      val project = projectRepo.insert(ProjectEntity.create("MemberTest", "did:admin"))
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      val added = projectRepo.addMember(project.id, biosample.id)
      assert(added)

      assert(projectRepo.isMember(project.id, biosample.id))
      assertEquals(projectRepo.countMembers(project.id), 1L)
    }
  }

  testTransactor.test("addMember is idempotent") { case (db, tx) =>
    tx.readWrite {
      val project = projectRepo.insert(ProjectEntity.create("IdempotentTest", "did:admin"))
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      // Add twice
      projectRepo.addMember(project.id, biosample.id)
      projectRepo.addMember(project.id, biosample.id)

      // Should still only have one member
      assertEquals(projectRepo.countMembers(project.id), 1L)
    }
  }

  testTransactor.test("removeMember removes biosample from project") { case (db, tx) =>
    tx.readWrite {
      val project = projectRepo.insert(ProjectEntity.create("RemoveTest", "did:admin"))
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      projectRepo.addMember(project.id, biosample.id)
      assert(projectRepo.isMember(project.id, biosample.id))

      val removed = projectRepo.removeMember(project.id, biosample.id)
      assert(removed)
      assert(!projectRepo.isMember(project.id, biosample.id))
    }
  }

  testTransactor.test("getMemberIds returns all member IDs") { case (db, tx) =>
    tx.readWrite {
      val project = projectRepo.insert(ProjectEntity.create("MembersTest", "did:admin"))
      val bs1 = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val bs2 = biosampleRepo.insert(BiosampleEntity.create("BS002", "D2"))
      val bs3 = biosampleRepo.insert(BiosampleEntity.create("BS003", "D3"))

      projectRepo.addMember(project.id, bs1.id)
      projectRepo.addMember(project.id, bs2.id)
      projectRepo.addMember(project.id, bs3.id)

      val memberIds = projectRepo.getMemberIds(project.id)
      assertEquals(memberIds.size, 3)
      assert(memberIds.contains(bs1.id))
      assert(memberIds.contains(bs2.id))
      assert(memberIds.contains(bs3.id))
    }
  }

  testTransactor.test("getMemberships returns memberships with timestamps") { case (db, tx) =>
    tx.readWrite {
      val project = projectRepo.insert(ProjectEntity.create("MembershipTest", "did:admin"))
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      projectRepo.addMember(project.id, biosample.id)

      val memberships = projectRepo.getMemberships(project.id)
      assertEquals(memberships.size, 1)
      assertEquals(memberships.head.projectId, project.id)
      assertEquals(memberships.head.biosampleId, biosample.id)
      assert(memberships.head.addedAt != null)
    }
  }

  testTransactor.test("getProjectsForBiosample returns projects containing biosample") { case (db, tx) =>
    tx.readWrite {
      val project1 = projectRepo.insert(ProjectEntity.create("Project1", "did:admin"))
      val project2 = projectRepo.insert(ProjectEntity.create("Project2", "did:admin"))
      val project3 = projectRepo.insert(ProjectEntity.create("Project3", "did:admin"))
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      projectRepo.addMember(project1.id, biosample.id)
      projectRepo.addMember(project2.id, biosample.id)
      // Not added to project3

      val projects = projectRepo.getProjectsForBiosample(biosample.id)
      assertEquals(projects.size, 2)
      assert(projects.contains(project1.id))
      assert(projects.contains(project2.id))
      assert(!projects.contains(project3.id))
    }
  }

  testTransactor.test("setMembers replaces all members") { case (db, tx) =>
    tx.readWrite {
      val project = projectRepo.insert(ProjectEntity.create("SetTest", "did:admin"))
      val bs1 = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val bs2 = biosampleRepo.insert(BiosampleEntity.create("BS002", "D2"))
      val bs3 = biosampleRepo.insert(BiosampleEntity.create("BS003", "D3"))

      // Add initial members
      projectRepo.addMember(project.id, bs1.id)
      projectRepo.addMember(project.id, bs2.id)
      assertEquals(projectRepo.countMembers(project.id), 2L)

      // Replace with new set
      projectRepo.setMembers(project.id, List(bs2.id, bs3.id))

      val memberIds = projectRepo.getMemberIds(project.id)
      assertEquals(memberIds.size, 2)
      assert(!memberIds.contains(bs1.id)) // Removed
      assert(memberIds.contains(bs2.id))  // Kept
      assert(memberIds.contains(bs3.id))  // Added
    }
  }

  testTransactor.test("addMember marks synced project as modified") { case (db, tx) =>
    tx.readWrite {
      val project = projectRepo.insert(ProjectEntity.create("SyncTest", "did:admin"))
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      // Mark project as synced
      projectRepo.markSynced(project.id, "at://test/1", "cid1")
      assertEquals(projectRepo.findById(project.id).get.meta.syncStatus, SyncStatus.Synced)

      // Add member - should mark as modified
      projectRepo.addMember(project.id, biosample.id)

      assertEquals(projectRepo.findById(project.id).get.meta.syncStatus, SyncStatus.Modified)
    }
  }

  testTransactor.test("unique constraint on project_name") { case (db, tx) =>
    val result = tx.readWrite {
      projectRepo.insert(ProjectEntity.create("DupeProject", "did:admin1"))
      projectRepo.insert(ProjectEntity.create("DupeProject", "did:admin2"))
    }

    assert(result.isLeft, "Should fail on duplicate project name")
  }
