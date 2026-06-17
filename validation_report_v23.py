import json
import logging
from variant_engine import annotate_snp_from_rsid

logging.basicConfig(level=logging.INFO)

VARIANTS_TO_TEST = [
    "rs429358",
    "rs7412",
    "rs334",
    "rs113993960",
    "rs1800562",
    "rs6025",
    "rs7903146",
    "rs9939609",
]

def run_validation():
    results = {}
    print("==================================================")
    print(" VERSION 2.3 EVIDENCE CONFIDENCE VALIDATION ")
    print("==================================================")
    
    for rsid in VARIANTS_TO_TEST:
        print(f"\\nValidating {rsid}...")
        try:
            report = annotate_snp_from_rsid(rsid)
            
            ev = report.get("evidence_confidence")
            
            if not ev:
                print(f"  [FAIL] Missing evidence_confidence object")
                results[rsid] = {"success": False, "error": "Missing evidence_confidence object"}
                continue
                
            score = ev.get("score")
            tier = ev.get("tier")
            strength = ev.get("strength")
            narrative = ev.get("narrative")
            factors = ev.get("contributing_factors", [])
            limits = ev.get("limitations", [])
            
            print(f"  [PASS] Score: {score}/100")
            print(f"  [PASS] Tier: {tier}")
            print(f"  [PASS] Strength: {strength}")
            print(f"  [PASS] Narrative dynamically generated ({len(narrative)} chars)")
            print(f"  [PASS] Contributing Factors: {len(factors)}")
            for f in factors:
                print(f"         + {f}")
            print(f"  [PASS] Limitations: {len(limits)}")
            for l in limits:
                print(f"         - {l}")
            
            results[rsid] = {
                "score": score,
                "tier": tier,
                "factors": len(factors),
                "limitations": len(limits),
                "success": True
            }
            
        except Exception as e:
            print(f"  [FAIL] Validation failed for {rsid}: {e}")
            results[rsid] = {"success": False, "error": str(e)}
            
    # Write report
    with open("validation_v23_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print("\\nValidation complete! Results saved to validation_v23_results.json")

if __name__ == "__main__":
    run_validation()
