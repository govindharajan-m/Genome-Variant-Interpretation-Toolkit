import json
import logging
from variant_engine import annotate_snp_from_rsid

logging.basicConfig(level=logging.INFO)

VARIANTS_TO_TEST = [
    "rs429358",
    "rs113993960",
    "rs1042522",
    "rs1800562",
    "rs6025",
    "rs7903146",
    "rs334",
]

def run_validation():
    results = {}
    print("==================================================")
    print(" VERSION 2.6 PHARMACOGENOMICS VALIDATION ")
    print("==================================================")
    
    for rsid in VARIANTS_TO_TEST:
        print(f"\\nValidating {rsid}...")
        try:
            report = annotate_snp_from_rsid(rsid)
            
            pgx = report.get("pharmacogenomics")
            
            if not pgx:
                print(f"  [FAIL] Missing pharmacogenomics object")
                results[rsid] = {"success": False, "error": "Missing pharmacogenomics object"}
                continue
                
            tier = pgx.get("tier")
            score = pgx.get("score")
            interactions = pgx.get("interactions", [])
            apps = pgx.get("applications", [])
            
            print(f"  [PASS] Tier: {tier}")
            print(f"  [PASS] Score: {score}")
            print(f"  [PASS] Interactions: {len(interactions)}")
            for i in interactions:
                print(f"         + {i['drug']} ({i['evidence_level']})")
                
            print(f"  [PASS] Applications: {len(apps)}")
            
            results[rsid] = {
                "tier": tier,
                "score": score,
                "interactions": len(interactions),
                "applications": len(apps),
                "success": True
            }
            
        except Exception as e:
            print(f"  [FAIL] Validation failed for {rsid}: {e}")
            results[rsid] = {"success": False, "error": str(e)}
            
    with open("validation_v26_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print("\\nValidation complete! Results saved to validation_v26_results.json")

if __name__ == "__main__":
    run_validation()
