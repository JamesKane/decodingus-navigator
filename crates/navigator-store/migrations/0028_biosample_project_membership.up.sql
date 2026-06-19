-- M:N Subject↔Project membership (FTDNA project-import design §4.1). A Subject (FTDNA kit) is the
-- canonical person and can belong to many projects; an admin runs many projects and the same kit
-- recurs across them. `biosample.project_id` is kept as a nullable "home project" pointer for one
-- release; `biosample_project` is the source of truth for membership going forward.
CREATE TABLE biosample_project (
    biosample_guid TEXT NOT NULL REFERENCES biosample(guid),
    project_id     INTEGER NOT NULL REFERENCES project(id),
    role           TEXT,          -- optional subgroup/branch label within the project
    added_at       TEXT NOT NULL, -- ISO-8601
    PRIMARY KEY (biosample_guid, project_id)
);
CREATE INDEX idx_biosample_project_proj ON biosample_project(project_id);

-- Backfill one membership row per existing single-project biosample.
INSERT INTO biosample_project (biosample_guid, project_id, role, added_at)
    SELECT guid, project_id, NULL, strftime('%Y-%m-%dT%H:%M:%SZ', 'now')
    FROM biosample
    WHERE project_id IS NOT NULL;
