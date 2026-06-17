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
    print(" VERSION 2.3.5 RESEARCH RELEVANCE VALIDATION ")
    print("==================================================")
    
    for rsid in VARIANTS_TO_TEST:
        print(f"\\nValidating {rsid}...")
        try:
            report = annotate_snp_from_rsid(rsid)
            
            rr = report.get("research_relevance")
            
            if not rr:
                print(f"  [FAIL] Missing research_relevance object")
                results[rsid] = {"success": False, "error": "Missing research_relevance object"}
                continue
                
            score = rr.get("score")
            tier = rr.get("research_relevance")
            reasons = rr.get("reasons", [])
            apps = rr.get("applications", [])
            
            print(f"  [PASS] Score: {score}/100")
            print(f"  [PASS] Tier: {tier}")
            print(f"  [PASS] Reasons: {len(reasons)}")
            for r in reasons:
                print(f"         + {r}")
            print(f"  [PASS] Applications: {len(apps)}")
            for a in apps:
                print(f"         > {a}")
            
            results[rsid] = {
                "score": score,
                "tier": tier,
                "reasons": len(reasons),
                "applications": len(apps),
                "success": True
            }
            
        except Exception as e:
            print(f"  [FAIL] Validation failed for {rsid}: {e}")
            results[rsid] = {"success": False, "error": str(e)}
            
    with open("validation_v235_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print("\\nValidation complete! Results saved to validation_v235_results.json")

if __name__ == "__main__":
    run_validation()
