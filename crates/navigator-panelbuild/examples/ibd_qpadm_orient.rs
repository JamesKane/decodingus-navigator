//! Do the IBD panel's chm13 REF/ALT and the qpAdm panel's REF/ALT agree at shared sites?
use navigator_analysis::ancestry::AncestryPanel;
use navigator_analysis::ibd_panel::IbdPanel;
use std::collections::HashMap;
fn main() -> anyhow::Result<()> {
    let ibd = IbdPanel::from_bytes(&std::fs::read(std::env::args().nth(1).unwrap())?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let qp = AncestryPanel::from_bytes(&std::fs::read(std::env::args().nth(2).unwrap())?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let m: HashMap<(String,i64),(char,char)> = ibd.sites.iter().map(|s| ((s.chm13.contig.clone(), s.chm13.position),(s.chm13.reference,s.chm13.alternate))).collect();
    let (mut ov, mut same, mut swap, mut other)=(0,0,0,0);
    for s in &qp.sites { if let Some(&(r,a))=m.get(&(s.contig.clone(),s.position)) { ov+=1;
        if (s.reference_allele,s.alternate_allele)==(r,a){same+=1} else if (s.reference_allele,s.alternate_allele)==(a,r){swap+=1} else {other+=1} } }
    println!("overlap {ov}: same {same}, swapped {swap}, other {other}");
    Ok(())
}
