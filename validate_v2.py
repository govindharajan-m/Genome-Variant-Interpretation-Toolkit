import variant_engine

rsids = ['rs429358', 'rs7412', 'rs334', 'rs1800562', 'rs1042522', 'rs1805007', 'rs113993960', 'rs6025', 'rs7903146', 'rs9939609']

for rsid in rsids:
    res = variant_engine.annotate_snp_from_rsid(rsid)
    if res.get('found'):
        gene_obj = res.get('gene')
        gene_name = gene_obj.get('gene') if gene_obj else 'None'
        impact_obj = res.get('impact')
        impact_level = impact_obj.get('impact_level') if impact_obj else 'None'
        clinvar_sig = res.get('clinvar', {}).get('clinical_significance') if res.get('clinvar') else 'None'
        print(f"[OK] {rsid}: SUCCESS (Gene: {gene_name}, Impact: {impact_level}, ClinVar: {clinvar_sig})")
    else:
        print(f"[FAIL] {rsid}: FAILED - {res.get('error')}")
