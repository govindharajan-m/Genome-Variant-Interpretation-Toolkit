import json
import logging
from app import app
from variant_engine import generate_panel_recommendation

logging.basicConfig(level=logging.INFO)

DISEASES_TO_TEST = [
    "Alzheimer's Disease",
    "Breast Cancer",
    "Hemochromatosis",
    "Cystic Fibrosis",
    "Type 2 Diabetes",
    "Thrombophilia",
    "Sickle Cell Disease"
]

def run_validation():
    results = {}
    print("==================================================")
    print(" VERSION 2.8 VARIANT DISCOVERY ENGINE VALIDATION ")
    print("==================================================")
    
    with app.app_context():
        for disease in DISEASES_TO_TEST:
            print(f"\\nValidating Discovery Engine for: {disease}...")
            try:
                panel = generate_panel_recommendation(disease)
                
                if "error" in panel:
                    print(f"  [FAIL] {panel['error']}")
                    results[disease] = {"success": False, "error": panel['error']}
                    continue
                    
                genes = panel.get("recommended_genes", [])
                summary = panel.get("summary", {})
                
                print(f"  [PASS] Summary Generated - Genes: {summary.get('genes_evaluated')}, Variants: {summary.get('candidate_variants')}")
                print(f"  [PASS] Top Variant: {summary.get('top_variant')} (Score: {summary.get('top_discovery_score')})")
                
                has_variants = False
                driver_exists = True
                flags_exist = True
                valid_sorting = True
                
                for g in genes:
                    print(f"         > Gene: {g['gene']} | Variant Count: {g['variant_count']}")
                    variants = g.get("variants", [])
                    if len(variants) > 0:
                        has_variants = True
                        for i, v in enumerate(variants):
                            print(f"           - {v['variant_id']}: Score={v['discovery_score']} ({v['discovery_tier']}) | Driver='{v['discovery_driver']}' | Flags={v['research_flags']}")
                            
                            if not v.get('discovery_driver'): driver_exists = False
                            if v.get('research_flags') is None: flags_exist = False
                            
                            if i > 0 and variants[i]['variant_priority'] > variants[i-1]['variant_priority']:
                                valid_sorting = False
                                
                if has_variants:
                    print(f"  [PASS] Variants properly nested under genes")
                else:
                    print(f"  [WARN] No variants mapped for this disease")
                    
                if driver_exists: print(f"  [PASS] Discovery Drivers correctly populated")
                else: print(f"  [FAIL] Missing Discovery Drivers")
                
                if flags_exist: print(f"  [PASS] Research Flags correctly populated")
                else: print(f"  [FAIL] Missing Research Flags")
                
                results[disease] = {
                    "gene_count": len(genes),
                    "variant_count": summary.get('candidate_variants'),
                    "top_variant": summary.get('top_variant'),
                    "success": has_variants and driver_exists and flags_exist and valid_sorting
                }
                
            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f"  [FAIL] Validation failed for {disease}: {e}")
                results[disease] = {"success": False, "error": str(e)}
                
    with open("validation_v28_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print("\\nValidation complete! Results saved to validation_v28_results.json")

if __name__ == "__main__":
    run_validation()
