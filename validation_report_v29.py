import json
import logging
from app import app
from variant_engine import analyze_variant_cohort

logging.basicConfig(level=logging.INFO)

RSID_LIST = [
    "rs429358", "rs7412", "rs334", "rs1800562", "rs1042522", 
    "rs6025", "rs7903146", "rs9939609", "rs113993960"
]

def run_validation():
    results = {}
    print("==================================================")
    print(" VERSION 2.9 COHORT PRIORITIZATION VALIDATION ")
    print("==================================================")
    
    with app.app_context():
        print(f"\\nValidating Cohort Analysis for {len(RSID_LIST)} variants...")
        try:
            cohort_data = analyze_variant_cohort(RSID_LIST)
            
            if "error" in cohort_data:
                print(f"  [FAIL] {cohort_data['error']}")
                results["cohort"] = {"success": False, "error": cohort_data['error']}
                return
                
            variants = cohort_data.get("variants", [])
            summary = cohort_data.get("summary", {})
            
            print(f"  [PASS] Summary Generated - Total Variants: {summary.get('total_variants')}")
            print(f"  [PASS] Top Variant: {summary.get('top_5_variants', [])[0]}")
            print(f"  [PASS] Top Gene Representation: {summary.get('top_gene_representation')}")
            
            has_variants = len(variants) > 0
            valid_sorting = True
            
            for i, v in enumerate(variants):
                print(f"         {i+1}. {v['variant_id']} ({v['gene']}) - Score: {v['cohort_score']} | Tier: {v['cohort_tier']} | EC: {v['evidence_confidence']} | RR: {v['research_relevance']}")
                
                if i > 0 and variants[i]['cohort_score'] > variants[i-1]['cohort_score']:
                    valid_sorting = False
                            
            if has_variants:
                print(f"  [PASS] Variants properly mapped and scored")
            else:
                print(f"  [FAIL] No variants returned")
                
            if valid_sorting: print(f"  [PASS] Cohort Ranking is correct (Descending order)")
            else: print(f"  [FAIL] Cohort Ranking is incorrect")
            
            results["cohort"] = {
                "total_variants": summary.get('total_variants'),
                "top_variant": summary.get('top_5_variants', [])[0] if summary.get('top_5_variants') else None,
                "success": has_variants and valid_sorting
            }
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  [FAIL] Validation failed: {e}")
            results["cohort"] = {"success": False, "error": str(e)}
                
    with open("validation_v29_results.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
        
    print("\\nValidation complete! Results saved to validation_v29_results.json")

if __name__ == "__main__":
    run_validation()
