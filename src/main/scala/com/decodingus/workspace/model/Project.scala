package com.decodingus.workspace.model

/**
 * A genealogy or research project that aggregates multiple biosamples.
 * This is a first-class record in the Atmosphere Lexicon (com.decodingus.atmosphere.project).
 *
 * @param atUri         The AT URI of this project record
 * @param meta          Record metadata for tracking changes and sync
 * @param projectName   Name of the project (e.g., 'Smith Surname Project')
 * @param description   Goals and scope of the research
 * @param administrator The DID of the researcher managing this project
 * @param memberRefs    AT URIs of biosample records associated with this project
 */
case class Project(
                    atUri: Option[String],
                    meta: RecordMeta,
                    projectName: String,
                    description: Option[String] = None,
                    administrator: String,
                    memberRefs: List[String] = List.empty
                  )
