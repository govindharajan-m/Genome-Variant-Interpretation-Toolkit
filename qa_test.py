import requests
import sys
import json

BASE_URL = "http://localhost:5000"

def test_endpoint(name, method, url, payload=None):
    print(f"Testing {name} ({method} {url})...", end=" ")
    try:
        if method == "POST":
            response = requests.post(BASE_URL + url, json=payload, timeout=5)
        elif method == "GET":
            response = requests.get(BASE_URL + url, params=payload, timeout=5)
        
        if response.status_code == 200:
            print("PASS")
            return True
        else:
            print(f"FAIL (Status: {response.status_code})")
            print("Response:", response.text[:200])
            return False
    except Exception as e:
        print(f"FAIL (Error: {str(e)})")
        return False

def run_tests():
    success = True
    
    # Test 1: Coordinate-based SNP
    success &= test_endpoint("Analyze SNP (coords)", "POST", "/api/analyze-snp", 
                             {"chromosome": "17", "position": 7676154, "ref": "C", "alt": "G"})
    
    # Test 2: rsID-based SNP
    success &= test_endpoint("Analyze SNP (rsID)", "POST", "/api/analyze-snp-rsid", 
                             {"rsid": "rs104894392"})
    
    # Test 3: CNV Analysis
    success &= test_endpoint("Analyze CNV", "POST", "/api/analyze-cnv", 
                             {"chromosome": "1", "start": 1000000, "end": 2000000, "cnv_type": "Deletion", "copy_number": 1})
    
    # Test 4: Batch Analysis
    success &= test_endpoint("Batch Analysis", "POST", "/api/batch", 
                             {"rsids": ["rs104894392", "rs113488022"]})
    
    # Test 5: Panel Designer
    print("Testing Panel Designer (POST /api/panel)...", end=" ")
    r = requests.post(BASE_URL + "/api/panel", data={"disease": "Alzheimer's Disease"}, timeout=15)
    if r.status_code == 200: print("PASS")
    else: 
        print(f"FAIL ({r.status_code})")
        success = False

    # Test 6: Cohort Analysis
    print("Testing Cohort Analysis (POST /api/cohort)...", end=" ")
    r = requests.post(BASE_URL + "/api/cohort", data={"rsids": "rs104894392,rs113488022"}, timeout=15)
    if r.status_code == 200: print("PASS")
    else: 
        print(f"FAIL ({r.status_code})")
        success = False

    # Test 7: Comparative Analysis
    print("Testing Comparative Analysis (POST /api/compare)...", end=" ")
    r = requests.post(BASE_URL + "/api/compare", data={"rsids": "rs104894392,rs113488022"}, timeout=15)
    if r.status_code == 200: print("PASS")
    else: 
        print(f"FAIL ({r.status_code})")
        success = False
                             
    # Test 8: HTML Page Loading
    pages = ["/", "/single-variant", "/cnv-analysis", "/batch-analysis", "/panel_designer", "/cohort_analysis", "/comparative_analysis"]
    for page in pages:
        res = requests.get(BASE_URL + page)
        if res.status_code == 200:
            print(f"Testing Page Load ({page})... PASS")
        else:
            print(f"Testing Page Load ({page})... FAIL ({res.status_code})")
            success = False

    if success:
        print("\nAll functional tests PASSED.")
        sys.exit(0)
    else:
        print("\nSome functional tests FAILED.")
        sys.exit(1)

if __name__ == "__main__":
    run_tests()
