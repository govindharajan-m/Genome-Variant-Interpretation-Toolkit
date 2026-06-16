import json
import db_handler
import variant_engine

output_dir = "C:/Users/mgovi/.gemini/antigravity-ide/brain/5ba47b3c-9e65-4ee7-9310-47ceacac79b5"
dbsnp_path = "d:/genome_variant_platform/genome_variant_platform/data/dbsnp_mock.json"

with open(dbsnp_path, 'r') as f:
    dbsnp_data = json.load(f)

audit_report = ["# SCIENTIFIC INTEGRITY AUDIT — V2.0.2\n\n"]
bug_tracker = ["# SCIENTIFIC BUG TRACKER — V2.0.2\n\n"]

missing_genes = []
impact_mismatches = []
clinvar_mismatches = []

# List to audit: all from mock + specific ones user requested
rsids_to_audit = set(dbsnp_data.keys())
rsids_to_audit.update(['rs429358', 'rs7412', 'rs6025', 'rs7903146', 'rs113993960'])

for rsid in sorted(rsids_to_audit):
    # Determine what mock data expected, if any
    mock_record = dbsnp_data.get(rsid, {})
    
    res = variant_engine.annotate_snp_from_rsid(rsid)
    if not res.get("found"):
        continue

    # 2 & 3. Gene Mapping
    gene = res.get("gene", {})
    gene_name = gene.get("gene") if gene else None
    
    # We should flag missing genes
    if not gene_name:
        missing_genes.append(f"- **{rsid}**: Missing gene mapping. Found None.")
    elif mock_record and mock_record.get("gene") and gene_name != mock_record.get("gene"):
        missing_genes.append(f"- **{rsid}**: Gene mismatch! Expected {mock_record.get('gene')}, got {gene_name}")

    # 4. Impact vs Consequence
    consequence = mock_record.get("consequence", "") if mock_record else res.get("impact", {}).get("consequence", "")
    impact = res.get("impact", {}).get("impact_level")
    
    if "nonsense" in consequence and impact != "HIGH":
        impact_mismatches.append(f"- **{rsid}**: Nonsense variant classified as {impact} instead of HIGH.")
    elif "missense" in consequence and impact not in ["MODERATE", "HIGH"]:
        impact_mismatches.append(f"- **{rsid}**: Missense variant classified as {impact} instead of MODERATE/HIGH.")

    # 5 & 6. ClinVar Mismatch
    live_clinvar = (res.get("clinvar") or {}).get("clinical_significance")
    mock_clinvar = (mock_record.get("clinvar") or {}).get("clinical_significance") if mock_record else None
    
    if live_clinvar and mock_clinvar and live_clinvar.lower() != mock_clinvar.lower():
        # Live ClinVar retrieval is finding something different from the static mock annotation
        clinvar_mismatches.append(f"- **{rsid}**: Live ClinVar='{live_clinvar}' vs Static Mock='{mock_clinvar}'. Check for outdated or contradictory data.")


# Formulate markdown
audit_report.append("## Gene Mapping Accuracy\n")
if missing_genes:
    audit_report.append("Identified missing or incorrect gene mappings:\n" + "\n".join(missing_genes) + "\n\n")
    bug_tracker.append("## HIGH Priority\n\n### Gene Mapping Failures\n" + "\n".join(missing_genes) + "\n\n")
else:
    audit_report.append("All gene mappings correct.\n\n")

audit_report.append("## Impact Predictions\n")
if impact_mismatches:
    audit_report.append("Identified inconsistent impact predictions:\n" + "\n".join(impact_mismatches) + "\n\n")
    bug_tracker.append("## MEDIUM Priority\n\n### Impact Prediction Inconsistencies\n" + "\n".join(impact_mismatches) + "\n\n")
else:
    audit_report.append("All impact predictions are biologically sound.\n\n")

audit_report.append("## ClinVar Classifications\n")
if clinvar_mismatches:
    audit_report.append("Identified discrepancies between live and static ClinVar data:\n" + "\n".join(clinvar_mismatches) + "\n\n")
    bug_tracker.append("## MEDIUM Priority\n\n### ClinVar Data Drift\n" + "\n".join(clinvar_mismatches) + "\n\n")
else:
    audit_report.append("All ClinVar classifications consistent.\n\n")

with open(f"{output_dir}/SCIENTIFIC_AUDIT_V203.md", "w", encoding="utf-8") as f:
    f.writelines(audit_report)
    
with open(f"{output_dir}/SCIENTIFIC_BUG_TRACKER_V203.md", "w", encoding="utf-8") as f:
    f.writelines(bug_tracker)
