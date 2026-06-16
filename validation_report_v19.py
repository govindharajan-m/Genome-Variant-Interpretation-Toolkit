import json
import json
from variant_engine import annotate_snp_from_rsid

def generate_validation_report():
    print("Generating validation report for Version 1.9 (GWAS Evidence Engine)...")
    report_lines = [
        "# Version 1.9 GWAS Evidence Validation Report",
        "",
        "This report verifies the extraction and automatic interpretation of GWAS associations.",
        ""
    ]

    # Required test cases
    test_cases = [
        {"rsid": "rs9939609", "desc": "FTO variant (Obesity/BMI)"},
        {"rsid": "rs429358", "desc": "APOE e4 variant (Alzheimer's)"},
        {"rsid": "rs7412", "desc": "APOE e2 variant"},
        {"rsid": "rs1042522", "desc": "TP53 variant"}
    ]

    for tc in test_cases:
        rsid = tc["rsid"]
        report_lines.append(f"## Testing {rsid} - {tc['desc']}")
        # Use variant engine to get processed interpretation
        record = annotate_snp_from_rsid(rsid)
        
        if record.get("error"):
            report_lines.append(f"Error: {record.get('error')}\n")
            continue

        gwas_assocs = record.get("gwas_associations")
        if not gwas_assocs:
            report_lines.append("**Result:** No GWAS associations available (Handled gracefully).\n")
            continue

        report_lines.append(f"**Total Associations Stored Internally:** {len(gwas_assocs)}")
        
        top_hit = record.get("gwas_top_hit")
        interp = record.get("gwas_interpretation")
        
        if top_hit:
            report_lines.append("\n### Top Association Displayed:")
            report_lines.append(f"- **Trait:** {top_hit.get('trait') or 'N/A'}")
            report_lines.append(f"- **P-value:** {top_hit.get('p_value')}")
            
            beta = top_hit.get('beta')
            or_val = top_hit.get('odds_ratio')
            if beta is not None:
                report_lines.append(f"- **Beta:** {beta}")
            if or_val is not None:
                report_lines.append(f"- **Odds Ratio:** {or_val}")
                
            report_lines.append(f"- **Study:** {top_hit.get('study') or 'N/A'}")
            report_lines.append(f"- **PMID:** {top_hit.get('pmid') or 'N/A'}")
            
        if interp:
            report_lines.append("\n### Interpretation Engine:")
            report_lines.append(f"- **Strength:** {interp.get('association_strength')}")
            report_lines.append(f"- **Narrative:** {interp.get('narrative')}")
            
        report_lines.append("\n---\n")

    with open("validation_report_v19.md", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))
    print("Validation report saved to validation_report_v19.md")

if __name__ == "__main__":
    generate_validation_report()
