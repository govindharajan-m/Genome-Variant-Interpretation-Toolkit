import json
from variant_engine import annotate_snp_from_rsid

VARIANTS = [
    "rs429358",
    "rs7412",
    "rs113993960"
]

def run():
    print("==================================================")
    print("V181 POPULATION REFINEMENT VALIDATION")
    print("==================================================")
    results = []
    
    for rsid in VARIANTS:
        print(f"\\nValidating {rsid}...")
        res = annotate_snp_from_rsid(rsid)
        
        comp = res.get("frequency_comparison")
        interp = res.get("frequency_interpretation")
        
        print(f"Interpretation: {interp}")
        if comp:
            print(f"SAS: {comp.get('south_asian')}")
            print(f"EUR: {comp.get('european')}")
            print(f"Abs Diff: {comp.get('absolute_difference')}")
            print(f"Fold Diff: {comp.get('fold_difference')}")
            print(f"Higher In: {comp.get('higher_in')}")
        else:
            print("No comparison data available.")
            
        results.append({
            "rsid": rsid,
            "comparison": comp
        })

    with open("validation_v181.json", "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    run()
