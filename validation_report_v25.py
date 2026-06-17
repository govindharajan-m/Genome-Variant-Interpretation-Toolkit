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
    print(" VERSION 2.5 ACMG ENGINE VALIDATION ")
    print("==================================================")
    
    for rsid in VARIANTS_TO_TEST:
        print(f"\\nValidating {rsid}...")
        try:
            report = annotate_snp_from_rsid(rsid)
            
            acmg = report.get("acmg_evidence")
            
            if not acmg:
                print(f"  [FAIL] Missing acmg_evidence object")
                results[rsid] = {"success": False, "error": "Missing acmg_evidence object"}
                continue
                
            cls = acmg.get("classification")
            trig = acmg.get("triggered_criteria", [])
            conf = acmg.get("mapping_confidence")
            
            print(f"  [PASS] Status: {cls}")
            print(f"  [PASS] Confidence: {conf}")
            print(f"  [PASS] Triggers: {', '.join(trig)}")
            
            for c in acmg.get("criteria", []):
                print(f"         + {c['code']} ({c['strength']}) - {c['reason']}")
            
            results[rsid] = {
                "status": cls,
                "confidence": conf,
                "triggers": trig,
                "success": True
            }
            
        except Exception as e:
            print(f"  [FAIL] Validation failed for {rsid}: {e}")
            results[rsid] = {"success": False, "error": str(e)}
            
    with open("validation_v25_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print("\\nValidation complete! Results saved to validation_v25_results.json")

if __name__ == "__main__":
    run_validation()
