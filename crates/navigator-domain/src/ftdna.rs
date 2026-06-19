//! FTDNA project-export parsers (FTDNA project-import design §3). Pure text→typed rows, no IO — the
//! app layer reads the files and hands the text here.
//!
//! Covers the two batch report CSVs that seed the importer's spine (Phase 1):
//! - `Member_Information` — the roster (§3.1)
//! - `Paternal_Ancestry` / `Maternal_Ancestry` — MDKA + clade path (§3.2, identical layout)
//!
//! The wide `YDNA_Results_Overview` Y-STR chart (§3.3) is parsed by [`crate::strprofile`]; this
//! module only handles the roster + ancestry files. All fields are looked up **by header name**
//! (not fixed position) since exports vary, columns are quoted-with-commas, and headers carry HTML
//! entities (`&gt;`, `&darr;`, `&amp;`) that are normalized here.

/// One member from `Member_Information` (§3.1). PII fields (`name`) are carried so the matcher can
/// fuzzy-compare, but they are never federated.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MemberRow {
    pub kit_number: String,
    pub name: Option<String>,
    /// `Access Granted` — pose-as gate + Big Y data tier (`Advanced`/`Limited`/…).
    pub access_granted: Option<String>,
    /// `Publicly Share DNA Results` (YES/NO) — federation consent.
    pub publicly_shares: Option<bool>,
}

/// One row from a `Paternal_Ancestry` / `Maternal_Ancestry` export (§3.2). The MDKA source.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct AncestryRow {
    pub kit_number: String,
    /// `Sub Group` — the project's clade/branch path (HTML-unescaped), e.g. `CTS4466>S1115>…`.
    pub sub_group: Option<String>,
    pub country: Option<String>,
    /// `Paternal/Maternal Ancestor Name` with the inline `b.`/`d.` dates stripped to [`Self::birth_year`]/
    /// [`Self::death_year`]; the leading name portion is kept here.
    pub ancestor_name: Option<String>,
    pub birth_year: Option<i32>,
    pub death_year: Option<i32>,
    /// `Map Location` (the `"No Location Saved"` sentinel dropped).
    pub origin_place: Option<String>,
    pub latitude: Option<f64>,
    pub longitude: Option<f64>,
}

/// Parse a `Member_Information` CSV. Skips rows whose `Kit Number` is blank.
pub fn parse_member_information(text: &str) -> Result<Vec<MemberRow>, String> {
    let (headers, mut rdr) = open(text)?;
    let kit = col(&headers, "Kit Number").ok_or("Member_Information: missing 'Kit Number' column")?;
    let name = col(&headers, "Name");
    let access = col(&headers, "Access Granted");
    let shares = col(&headers, "Publicly Share DNA Results");

    let mut out = Vec::new();
    for rec in rdr.records() {
        let rec = rec.map_err(|e| e.to_string())?;
        let kit_number = get(&rec, kit);
        if !is_real_kit(&kit_number) {
            continue;
        }
        out.push(MemberRow {
            kit_number,
            name: name.and_then(|i| nonblank(get(&rec, i))),
            access_granted: access.and_then(|i| nonblank(get(&rec, i))),
            publicly_shares: shares.and_then(|i| parse_yes_no(&get(&rec, i))),
        });
    }
    Ok(out)
}

/// Parse a `Paternal_Ancestry` / `Maternal_Ancestry` CSV (same layout). Skips blank-kit rows.
pub fn parse_ancestry(text: &str) -> Result<Vec<AncestryRow>, String> {
    let (headers, mut rdr) = open(text)?;
    let kit = col(&headers, "Kit Number").ok_or("Ancestry: missing 'Kit Number' column")?;
    let sub_group = col(&headers, "Sub Group");
    let country = col(&headers, "Country");
    // "Paternal Ancestor Name" or "Maternal Ancestor Name".
    let ancestor = headers.iter().position(|h| h.ends_with("Ancestor Name"));
    let map_loc = col(&headers, "Map Location");
    let lat = col(&headers, "Latitude");
    let lon = col(&headers, "Longitude");

    let mut out = Vec::new();
    for rec in rdr.records() {
        let rec = rec.map_err(|e| e.to_string())?;
        let kit_number = get(&rec, kit);
        if !is_real_kit(&kit_number) {
            continue;
        }
        let (name, birth_year, death_year) = match ancestor.and_then(|i| nonblank(get(&rec, i))) {
            Some(raw) => parse_ancestor_name(&raw),
            None => (None, None, None),
        };
        out.push(AncestryRow {
            kit_number,
            sub_group: sub_group.and_then(|i| nonblank_clade(get(&rec, i))),
            country: country
                .and_then(|i| nonblank(get(&rec, i)))
                .filter(|c| c != "Unknown Origin"),
            ancestor_name: name,
            birth_year,
            death_year,
            origin_place: map_loc
                .and_then(|i| nonblank(get(&rec, i)))
                .filter(|p| p != "No Location Saved"),
            latitude: lat.and_then(|i| parse_coord(&get(&rec, i))),
            longitude: lon.and_then(|i| parse_coord(&get(&rec, i))),
        });
    }
    Ok(out)
}

// ---- helpers --------------------------------------------------------------

/// Build a `csv::Reader` over the text with the header row pulled out + HTML-unescaped. Flexible
/// (some FTDNA exports have ragged trailing columns) and trims so quoted, space-padded cells clean up.
fn open(text: &str) -> Result<(Vec<String>, csv::Reader<&[u8]>), String> {
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .has_headers(true)
        .trim(csv::Trim::All)
        .from_reader(text.as_bytes());
    let headers = rdr
        .headers()
        .map_err(|e| e.to_string())?
        .iter()
        .map(unescape_html)
        .collect();
    Ok((headers, rdr))
}

/// First column index whose (unescaped) header equals `name`.
fn col(headers: &[String], name: &str) -> Option<usize> {
    headers.iter().position(|h| h == name)
}

fn get(rec: &csv::StringRecord, idx: usize) -> String {
    unescape_html(rec.get(idx).unwrap_or("").trim())
}

fn nonblank(s: String) -> Option<String> {
    let t = s.trim();
    if t.is_empty() || t == "-" {
        None
    } else {
        Some(t.to_string())
    }
}

/// A `Sub Group` value is only a clade path when it actually contains a lineage (`>`); FTDNA also
/// uses free-text placeholders there ("Not Yet Tested Positive for Relevant SNPs").
fn nonblank_clade(s: String) -> Option<String> {
    nonblank(s).filter(|v| v.contains('>'))
}

/// A real kit number is non-empty and not one of the two leading non-member sentinel rows
/// (`00000.` panel / `MIN`) the Y-STR overview carries — harmless to guard here too.
fn is_real_kit(kit: &str) -> bool {
    let k = kit.trim();
    !k.is_empty() && k != "MIN" && !k.starts_with("00000")
}

/// YES/NO (case-insensitive) → bool; anything else → None.
fn parse_yes_no(s: &str) -> Option<bool> {
    match s.trim().to_ascii_uppercase().as_str() {
        "YES" => Some(true),
        "NO" => Some(false),
        _ => None,
    }
}

/// Coordinate cell → f64, dropping the FTDNA `0` sentinel (means "no location").
fn parse_coord(s: &str) -> Option<f64> {
    let v: f64 = s.trim().parse().ok()?;
    (v != 0.0).then_some(v)
}

/// Minimal HTML-entity unescape for the entities FTDNA emits in headers/values.
fn unescape_html(s: &str) -> String {
    s.replace("&gt;", ">")
        .replace("&lt;", "<")
        .replace("&amp;", "&")
        .replace("&darr;", "")
        .replace("&#39;", "'")
        .replace("&quot;", "\"")
        .trim()
        .to_string()
}

/// Split an FTDNA ancestor field into `(name, birth_year, death_year)`. The dates are embedded
/// inline in varied shapes — `"Thomas Michael Kane, b. 1830 Clare, IE d. 1908 WI"`,
/// `"Joseph Abbett, b. 19 Mar 1819 and d. 2 Nov 1852"` — so we locate `b.`/`d.` markers and take the
/// first 4-digit year after each. The name is everything before the first marker (trailing comma
/// trimmed).
fn parse_ancestor_name(raw: &str) -> (Option<String>, Option<i32>, Option<i32>) {
    let lower = raw.to_ascii_lowercase();
    let b_pos = find_marker(&lower, "b.");
    let d_pos = find_marker(&lower, "d.");

    let birth = b_pos.and_then(|p| first_year(&raw[p..d_pos.unwrap_or(raw.len()).max(p)]));
    let death = d_pos.and_then(|p| first_year(&raw[p..]));

    let name_end = [b_pos, d_pos].iter().flatten().copied().min().unwrap_or(raw.len());
    let name = raw[..name_end].trim().trim_end_matches(',').trim();
    let name = if name.is_empty() { None } else { Some(name.to_string()) };
    (name, birth, death)
}

/// Byte offset of a `b.`/`d.` date marker, requiring a word boundary before it (so the `b` in
/// "Abbett" doesn't match). Returns the offset of the marker letter.
fn find_marker(lower: &str, marker: &str) -> Option<usize> {
    let bytes = lower.as_bytes();
    let mut from = 0;
    while let Some(rel) = lower[from..].find(marker) {
        let at = from + rel;
        let boundary = at == 0 || !bytes[at - 1].is_ascii_alphanumeric();
        if boundary {
            return Some(at);
        }
        from = at + 1;
    }
    None
}

/// First 4-digit run in `s` read as a plausible year (1000–2999).
fn first_year(s: &str) -> Option<i32> {
    let bytes = s.as_bytes();
    let mut i = 0;
    while i + 4 <= bytes.len() {
        if bytes[i].is_ascii_digit() {
            let run_len = bytes[i..].iter().take_while(|b| b.is_ascii_digit()).count();
            if run_len == 4 {
                if let Ok(y) = s[i..i + 4].parse::<i32>() {
                    if (1000..3000).contains(&y) {
                        return Some(y);
                    }
                }
            }
            i += run_len;
        } else {
            i += 1;
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn member_information_parses_kit_access_consent() {
        let csv = "Kit Number,Family Tree,Name,Email,Note,Release,Kit Back,Last Sign In,Access Granted,Allows MyHeritage Connection,Publicly Share DNA Results,Remove From Group\n\
                   B5163,YES,Jane Roe,jane@x.com,,YES,7/3/2013,6/3/2026,Limited,NO,YES,\n\
                   B92390,NO,John Doe,,note text,YES,N/A,7/22/2023,Full,NO,NO,\n\
                   ,NO,blank kit row,,,,,,,,,\n";
        let rows = parse_member_information(csv).unwrap();
        assert_eq!(rows.len(), 2, "blank-kit row skipped");
        assert_eq!(rows[0].kit_number, "B5163");
        assert_eq!(rows[0].access_granted.as_deref(), Some("Limited"));
        assert_eq!(rows[0].publicly_shares, Some(true));
        assert_eq!(rows[1].publicly_shares, Some(false));
        assert_eq!(rows[1].name.as_deref(), Some("John Doe"));
    }

    #[test]
    fn ancestry_parses_clade_dates_coords_with_quoting() {
        // B5163's real (redacted) shape: quoted ancestor with inline b./d., real coords, &gt; clade.
        let csv = "Kit Number,Name&darr;,Sub Group,Email,Country,Comment,Paternal Ancestor Name,Map Location,Latitude,Longitude,Family Tree,Family Tree,Remove From Group\n\
                   B5163,REDACTED,31050. CTS4466&gt;S1115&gt;FGC29071,,Ireland,,\"Thomas Michael Kane, b. 1830 Clare, IE d. 1908 WI\",\"Creegh South, Co. Clare, Ireland\",52.75,-9.43,,WikiTree,\n\
                   B625697,REDACTED,Not Yet Tested Positive for Relevant SNPs,,Unknown Origin,,James Joseph Dinn,No Location Saved,0,0,,WikiTree,\n";
        let rows = parse_ancestry(csv).unwrap();
        assert_eq!(rows.len(), 2);

        let b = &rows[0];
        assert_eq!(b.kit_number, "B5163");
        assert_eq!(b.sub_group.as_deref(), Some("31050. CTS4466>S1115>FGC29071"));
        assert_eq!(b.country.as_deref(), Some("Ireland"));
        assert_eq!(b.ancestor_name.as_deref(), Some("Thomas Michael Kane"));
        assert_eq!(b.birth_year, Some(1830));
        assert_eq!(b.death_year, Some(1908));
        assert_eq!(b.origin_place.as_deref(), Some("Creegh South, Co. Clare, Ireland"));
        assert_eq!(b.latitude, Some(52.75));
        assert_eq!(b.longitude, Some(-9.43));

        // Sentinels dropped; free-text Sub Group is not a clade.
        let o = &rows[1];
        assert_eq!(o.sub_group, None);
        assert_eq!(o.country, None, "Unknown Origin dropped");
        assert_eq!(o.origin_place, None, "No Location Saved dropped");
        assert_eq!(o.latitude, None, "0 sentinel dropped");
        assert_eq!(o.ancestor_name.as_deref(), Some("James Joseph Dinn"));
        assert_eq!(o.birth_year, None);
    }

    #[test]
    fn ancestor_dates_handle_day_month_year_form_and_word_boundary() {
        // "Abbett" must not trip the b. marker; "b. 19 Mar 1819 and d. 2 Nov 1852" → 1819/1852.
        let (name, b, d) = parse_ancestor_name("Joseph Abbett, b. 19 Mar 1819 and d. 2 Nov 1852");
        assert_eq!(name.as_deref(), Some("Joseph Abbett"));
        assert_eq!(b, Some(1819));
        assert_eq!(d, Some(1852));

        // No dates at all → whole string is the name.
        let (name, b, d) = parse_ancestor_name("James Joseph Dinn");
        assert_eq!(name.as_deref(), Some("James Joseph Dinn"));
        assert_eq!((b, d), (None, None));

        // Birth only.
        let (_, b, d) = parse_ancestor_name("Some One b. 1900");
        assert_eq!((b, d), (Some(1900), None));
    }
}
