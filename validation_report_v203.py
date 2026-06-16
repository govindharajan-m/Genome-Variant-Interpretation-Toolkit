import variant_engine

rsids = ['rs6025', 'rs7903146', 'rs429358', 'rs7412', 'rs334', 'rs1800562']

print("=== V2.0.3 Validation Report ===")
for rsid in rsids:
    res = variant_engine.annotate_snp_from_rsid(rsid)
    if not res.get("found"):
        print(f"[FAIL] {rsid}: NOT FOUND - {res.get('error')}")
        continue
    
    gene_obj = res.get("gene")
    gene = gene_obj.get("gene") if gene_obj else None
    
    impact = res.get("impact", {}).get("impact_level")
    clinvar = res.get("clinvar", {}).get("clinical_significance")
    
    # Just printing the summary of verification
    status = "SUCCESS" if gene else "FAILED (Missing Gene)"
    print(f"[{status}] {rsid}: Gene={gene}, Impact={impact}, ClinVar='{clinvar}'")

print("Validation complete.")
