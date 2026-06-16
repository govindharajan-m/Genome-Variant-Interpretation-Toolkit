import json
from variant_engine import annotate_snp_from_rsid

VARIANTS_TO_TEST = [
    "rs429358",
    "rs7412",
    "rs334",
    "rs1800562",
    "rs6025",
    "rs7903146",
    "rs9939609",
    "rs113993960"
]

def run_validation():
    print("==================================================")
    print("VERSION 2.1 PUBMED VALIDATION")
    print("==================================================")
    
    results = []
    
    for rsid in VARIANTS_TO_TEST:
        print(f"\\nValidating {rsid}...")
        res = annotate_snp_from_rsid(rsid)
        
        if not res.get("found"):
            print(f"  [ERROR] Not found in mock database!")
            continue
            
        gene = res.get("gene", {}).get("gene", "Unknown")
        pubmed = res.get("pubmed", {})
        
        print(f"  Gene: {gene}")
        print(f"  PubMed Available: {pubmed.get('available')}")
        print(f"  Paper Count: {pubmed.get('paper_count')}")
        print(f"  Evidence Level: {pubmed.get('evidence_level')}")
        print(f"  Top Papers Retrieved: {len(pubmed.get('papers', []))}")
        
        if pubmed.get('papers'):
            print(f"  Top PMID: {pubmed.get('papers')[0].get('pmid')}")
            
        results.append({
            "rsid": rsid,
            "gene": gene,
            "paper_count": pubmed.get("paper_count", 0),
            "evidence_level": pubmed.get("evidence_level", "No Evidence")
        })

    with open("validation_v21_results.json", "w") as f:
        json.dump(results, f, indent=2)
        
    print("\\n==================================================")
    print("VALIDATION COMPLETE")
    print("==================================================")

if __name__ == "__main__":
    run_validation()
