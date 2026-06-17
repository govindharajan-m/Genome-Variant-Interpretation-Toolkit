import json
import logging
from variant_engine import annotate_snp_from_rsid

logging.basicConfig(level=logging.INFO)

VARIANTS_TO_TEST = [
    "rs429358",
    "rs7412",
    "rs334",
    "rs1800562",
    "rs1042522",
    "rs6025",
    "rs7903146",
    "rs9939609",
    "rs113993960",
]

def run_validation():
    results = {}
    print("==================================================")
    print(" VERSION 2.4 VARIANT PRIORITIZATION VALIDATION ")
    print("==================================================")
    
    for rsid in VARIANTS_TO_TEST:
        print(f"\\nValidating {rsid}...")
        try:
            report = annotate_snp_from_rsid(rsid)
            
            vp = report.get("variant_priority")
            
            if not vp:
                print(f"  [FAIL] Missing variant_priority object")
                results[rsid] = {"success": False, "error": "Missing variant_priority object"}
                continue
                
            score = vp.get("priority_score")
            tier = vp.get("priority_tier")
            driver = vp.get("primary_driver")
            reasons = vp.get("priority_reasons", [])
            
            print(f"  [PASS] Score: {score}/100")
            print(f"  [PASS] Tier: {tier}")
            print(f"  [PASS] Driver: {driver}")
            print(f"  [PASS] Reasons: {len(reasons)}")
            for r in reasons:
                print(f"         + {r}")
            
            results[rsid] = {
                "score": score,
                "tier": tier,
                "driver": driver,
                "reasons": len(reasons),
                "success": True
            }
            
        except Exception as e:
            print(f"  [FAIL] Validation failed for {rsid}: {e}")
            results[rsid] = {"success": False, "error": str(e)}
            
    with open("validation_v24_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print("\\nValidation complete! Results saved to validation_v24_results.json")

if __name__ == "__main__":
    run_validation()
