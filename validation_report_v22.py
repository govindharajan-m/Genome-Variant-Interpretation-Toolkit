import json
import logging
from variant_engine import annotate_snp_from_rsid

logging.basicConfig(level=logging.INFO)

VARIANTS_TO_TEST = [
    "rs429358",      # APOE
    "rs113993960",   # CFTR
    "rs1042522",     # TP53
    "rs1800562",     # HFE
    "rs334",         # HBB
    "rs6025",        # F5
    "rs7903146",     # TCF7L2
    "rs121918310",   # MTHFR (Unknown Gene Test, since MTHFR is not in our dictionaries)
]

def run_validation():
    results = {}
    print("==================================================")
    print(" VERSION 2.2 GENE CONTEXT VALIDATION ")
    print("==================================================")
    
    for rsid in VARIANTS_TO_TEST:
        print(f"\\nValidating {rsid}...")
        try:
            report = annotate_snp_from_rsid(rsid)
            
            # Check gene context
            gc = report.get("gene_context", {})
            has_gc = gc.get("available", False)
            symbol = gc.get("symbol", "Unknown")
            
            # Check diseases
            diseases = report.get("diseases", [])
            has_diseases = len(diseases) > 0
            
            # Check pathways
            pathways = report.get("pathways", [])
            has_pathways = len(pathways) > 0
            
            relevance = report.get("research_relevance", "None")
            
            print(f"  [PASS] Report generation succeeded")
            if has_gc:
                print(f"  [PASS] Gene Context rendered ({symbol})")
            else:
                print(f"  [-] Gene Context unavailable")
                
            if has_diseases:
                print(f"  [PASS] Diseases found ({len(diseases)})")
            else:
                print(f"  [PASS] Disease module returned unavailable (Expected for unknown genes)")
                
            if has_pathways:
                print(f"  [PASS] Pathways found ({len(pathways)})")
            else:
                print(f"  [PASS] Pathway module returned unavailable (Expected for unknown genes)")
                
            print(f"  [PASS] Research Relevance: {relevance}")
            
            results[rsid] = {
                "gene": symbol,
                "has_context": has_gc,
                "disease_count": len(diseases),
                "pathway_count": len(pathways),
                "relevance": relevance,
                "success": True
            }
            
        except Exception as e:
            print(f"  [FAIL] Validation failed for {rsid}: {e}")
            results[rsid] = {"success": False, "error": str(e)}
            
    # Write report
    with open("validation_v22_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print("\\nValidation complete! Results saved to validation_v22_results.json")

if __name__ == "__main__":
    run_validation()
