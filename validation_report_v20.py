import json
from db_handler import get_dbsnp_record
from variant_engine import annotate_snp_from_rsid

def validate():
    test_rsids = ["rs429358", "rs7412", "rs334", "rs1805007", "rs6025", "rs9939609"]
    report = ["# Version 2.0 — Live ClinVar Validation Report\n"]
    
    for rsid in test_rsids:
        print(f"Validating {rsid}...")
        report.append(f"## Variant: {rsid}")
        
        # Test 1: Full payload extraction
        result = annotate_snp_from_rsid(rsid)
        if "error" in result:
            report.append(f"**Status:** ERROR - {result['error']}\n")
            continue
            
        cv = result.get("clinvar")
        if not cv:
            report.append("**ClinVar:** No annotation available.\n")
            continue
            
        report.append("### Extracted Data")
        report.append(f"- **Clinical Significance:** {cv.get('clinical_significance')}")
        report.append(f"- **Review Status:** {cv.get('review_status')}")
        report.append(f"- **Confidence Level:** {cv.get('confidence_level')}")
        report.append(f"- **Consensus Score:** {cv.get('consensus_score')}%")
        report.append(f"- **Accession:** {cv.get('accession')}")
        
        if cv.get('conflicting_interpretations'):
            report.append("\n### Conflicting Interpretations Found")
            for cls, count in cv.get('conflict_breakdown', {}).items():
                report.append(f"- {cls}: {count}")
                
        report.append("\n### Interpretation Summary Engine")
        report.append(f"> {result.get('interpretation_summary', 'N/A')}\n")
        
        report.append("---\n")
        
    with open("validation_report_v20.md", "w") as f:
        f.write("\n".join(report))
        
    print("Validation report generated: validation_report_v20.md")

if __name__ == "__main__":
    validate()
