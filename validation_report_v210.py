import requests
import json
import sys

def test_v210():
    print("Testing Comparative Variant Analysis Engine...")
    
    url = "http://localhost:5000/api/compare"
    test_rsids = "rs429358, rs7412, rs334, rs1800562, rs1042522, rs6025, rs7903146, rs9939609, rs113993960"
    
    try:
        r = requests.post(url, data={"rsids": test_rsids})
        if r.status_code != 200:
            print(f"[FAIL] Server returned {r.status_code}: {r.text}")
            sys.exit(1)
            
        data = r.json()
        
        # 1. Structure Verification
        if "variants" not in data or "summary" not in data:
            print("[FAIL] Missing variants or summary block in response")
            sys.exit(1)
            
        print("[PASS] Payload structure correct.")
            
        # 2. Variants Array Verification
        if len(data["variants"]) == 0:
            print("[FAIL] No variants evaluated")
            sys.exit(1)
            
        v = data["variants"][0]
        required_keys = ["variant_id", "gene", "evidence_confidence", "research_relevance", "variant_priority", "discovery_score", "cohort_score", "clin_score", "strength_matrix"]
        missing = [k for k in required_keys if k not in v]
        if missing:
            print(f"[FAIL] Variants missing required keys: {missing}")
            sys.exit(1)
            
        print(f"[PASS] Variant schema validated ({len(data['variants'])} variants processed).")
        
        # 3. Insights Verification
        summary = data["summary"]
        if not summary.get("insights"):
            print("[FAIL] No comparative insights generated")
            sys.exit(1)
            
        print("[PASS] Comparative insights generated successfully:")
        for ins in summary["insights"]:
            print(f"  -> {ins}")
            
        print("\nAll Tests Passed!")
        
    except Exception as e:
        print(f"Exception during testing: {e}")

if __name__ == "__main__":
    test_v210()
