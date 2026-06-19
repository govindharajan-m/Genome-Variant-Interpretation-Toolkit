import requests

def test_v296():
    print("Testing Security Enhancements...")
    
    # 1. Test Rate Limiting
    print("Testing Rate Limiting (/api/cohort)...")
    url = "http://localhost:5000/api/cohort"
    
    responses = []
    for _ in range(12):
        try:
            r = requests.post(url, data={})
            responses.append(r.status_code)
        except Exception:
            pass
            
    print(f"Status codes received: {responses}")
    if 429 in responses:
        print("[PASS] Rate limiting is ACTIVE.")
    else:
        print("[FAIL] Rate limiting not hit.")
        
    # 2. Test CSP Headers
    print("Testing CSP Headers...")
    try:
        r = requests.get("http://localhost:5000/")
        headers = r.headers
        if "Content-Security-Policy" in headers:
            print("[PASS] CSP Header found.")
        else:
            print("[FAIL] CSP Header missing.")
            
        if "X-Frame-Options" in headers and headers["X-Frame-Options"] == "DENY":
            print("[PASS] X-Frame-Options found.")
        else:
            print("[FAIL] X-Frame-Options missing.")
    except Exception:
        print("Could not connect to server to test headers.")

if __name__ == "__main__":
    test_v296()
