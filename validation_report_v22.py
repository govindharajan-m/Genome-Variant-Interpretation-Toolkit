import json
import logging
from variant_engine import annotate_snp_from_rsid
from db_handler import _fetch_live_gene_context, _fetch_gene_diseases, _fetch_gene_pathways

logging.basicConfig(level=logging.ERROR)

GENES_TO_TEST = ["APOE", "CFTR", "HBB", "HFE", "TP53", "F5", "TCF7L2"]

VARIANTS_TO_TEST = {
    "APOE": "rs429358",
    "CFTR": "rs113993960",
    "HBB": "rs334",
    "F5": "rs6025",
    "TCF7L2": "rs7903146"
}

def validate():
    print("="*60)
    print("VERSION 2.2 — GENE CONTEXT VALIDATION")
    print("="*60)
    
    print("\\n[Phase 1-3] Validating Live Retrieval Layers...")
    for gene in GENES_TO_TEST:
        print(f"\\nFetching Data for: {gene}")
        
        ctx = _fetch_live_gene_context(gene)
        dis = _fetch_gene_diseases(gene)
        pth = _fetch_gene_pathways(gene)
        
        print(f"  Context:   {'✅ SUCCESS' if ctx.get('available') else '❌ FAILED'} -> Name: {ctx.get('full_name')} | Source: {ctx.get('source')}")
        print(f"  Diseases:  {'✅ SUCCESS' if dis.get('available') else '❌ FAILED'} -> Found {len(dis.get('diseases', []))}")
        print(f"  Pathways:  {'✅ SUCCESS' if pth.get('available') else '❌ FAILED'} -> Found {len(pth.get('pathways', []))}")
        
    print("\\n" + "="*60)
    print("[Phase 4-5] Validating Engine Integration & Relevance Score...")
    print("="*60)
    
    for gene, rsid in VARIANTS_TO_TEST.items():
        print(f"\\nAnnotating {rsid} ({gene})...")
        res = annotate_snp_from_rsid(rsid)
        if res.get("found"):
            rel = res.get("research_relevance", "Unknown")
            has_ctx = res.get("gene_context", {}).get("available", False)
            dis_count = len(res.get("diseases", {}).get("diseases", []))
            
            print(f"  rsID: {rsid}")
            print(f"  Context Embedded: {'✅ YES' if has_ctx else '❌ NO'}")
            print(f"  Disease Count: {dis_count}")
            print(f"  Research Relevance: {rel}")
        else:
            print(f"  ❌ FAILED to annotate {rsid}")

if __name__ == "__main__":
    validate()
