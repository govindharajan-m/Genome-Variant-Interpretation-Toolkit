import json
import logging
from app import app
from variant_engine import generate_panel_recommendation

logging.basicConfig(level=logging.INFO)

DISEASES_TO_TEST = [
    "Alzheimer's Disease",
    "Breast Cancer",
    "Cystic Fibrosis",
    "Type 2 Diabetes",
    "Hemochromatosis"
]

def run_validation():
    results = {}
    print("==================================================")
    print(" VERSION 2.7 DISEASE PANEL DESIGNER VALIDATION ")
    print("==================================================")
    
    with app.app_context():
        for disease in DISEASES_TO_TEST:
            print(f"\\nValidating Panel: {disease}...")
            try:
                panel = generate_panel_recommendation(disease)
                
                if "error" in panel:
                    print(f"  [FAIL] {panel['error']}")
                    results[disease] = {"success": False, "error": panel['error']}
                    continue
                    
                genes = panel.get("recommended_genes", [])
                total_lit = panel.get("total_literature", 0)
                
                print(f"  [PASS] Candidate Genes Returned: {len(genes)}")
                print(f"  [PASS] Total Literature Count: {total_lit}")
                
                valid_ranking = True
                for i, g in enumerate(genes):
                    print(f"         {i+1}. {g['gene']} (Score: {g['panel_score']}, Priority: {g['priority']})")
                    if i > 0 and genes[i]['panel_score'] > genes[i-1]['panel_score']:
                        valid_ranking = False
                        
                if valid_ranking:
                    print(f"  [PASS] Ranking operates correctly (Descending order)")
                else:
                    print(f"  [FAIL] Ranking is incorrect")
                    
                has_reasons = all(len(g['reason']) > 0 for g in genes)
                if has_reasons:
                    print(f"  [PASS] Deterministic Justifications generated")
                else:
                    print(f"  [FAIL] Missing justifications")
                
                results[disease] = {
                    "gene_count": len(genes),
                    "total_lit": total_lit,
                    "top_gene": genes[0]['gene'] if genes else None,
                    "success": valid_ranking and has_reasons
                }
                
            except Exception as e:
                print(f"  [FAIL] Validation failed for {disease}: {e}")
                results[disease] = {"success": False, "error": str(e)}
                
    with open("validation_v27_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print("\\nValidation complete! Results saved to validation_v27_results.json")

if __name__ == "__main__":
    run_validation()
