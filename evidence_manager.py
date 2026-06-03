"""
evidence_manager.py — Evidence Retrieval & Scoring System
===========================================================
Unified framework to aggregate and score variant evidence from
multiple resources (dbSNP, ClinVar, GeneCards, PubMed).
"""

from db_handler import (
    get_dbsnp_record,
    get_clinvar_record,
    get_pubmed_mock,
    get_genecards_mock
)

def calculate_evidence_score(evidence: dict) -> tuple[int, str]:
    """
    Calculate an evidence confidence score based on available data.
    
    Returns:
        tuple of (score, strength_label)
    """
    score = 0
    
    if evidence.get("dbsnp"):
        score += 20
        
    if evidence.get("clinvar"):
        score += 30
        
    if evidence.get("genecards"):
        score += 20
        
    pubmed_count = len(evidence.get("pubmed", []))
    score += min(30, pubmed_count * 5)
    
    if score >= 70:
        strength = "High"
    elif score >= 40:
        strength = "Moderate"
    else:
        strength = "Limited"
        
    return score, strength


def build_evidence_object(rsid: str, gene: str) -> dict:
    """
    Construct a unified evidence object for a given variant.
    """
    evidence = {
        "variant": {
            "rsid": rsid,
            "gene": gene
        },
        "dbsnp": {},
        "clinvar": {},
        "genecards": {},
        "pubmed": [],
        "warnings": [],
        "evidence_score": 0,
        "evidence_strength": "Limited"
    }

    if not rsid:
        return evidence
        
    key = rsid.lower().strip()
    if not key.startswith("rs"):
        key = "rs" + key

    # 1. dbSNP
    dbsnp_record = get_dbsnp_record(key)
    if dbsnp_record:
        evidence["dbsnp"] = {
            "id": key,
            "url": f"https://www.ncbi.nlm.nih.gov/snp/{key}"
        }
    else:
        evidence["warnings"].append("Variant not found in dbSNP.")

    # 2. ClinVar
    clinvar_record = get_clinvar_record(key)
    if clinvar_record and clinvar_record.get("accession"):
        evidence["clinvar"] = {
            "id": f"ClinVar ID {clinvar_record['accession']}",
            "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_record['accession']}/",
            "classification": clinvar_record.get("clinical_significance", "Unknown"),
            "review_status": clinvar_record.get("review_status", "Unreviewed"),
            "disease_associations": [clinvar_record.get("condition")] if clinvar_record.get("condition") else []
        }
    else:
        evidence["warnings"].append("No ClinVar annotation available.")

    # 3. GeneCards
    if gene and gene != "—":
        genecards_info = get_genecards_mock(gene)
        if genecards_info:
            evidence["genecards"] = genecards_info
        else:
             evidence["warnings"].append(f"No GeneCards entry found for {gene}.")

    # 4. PubMed
    pubmed_refs = get_pubmed_mock(key, gene)
    if pubmed_refs:
        evidence["pubmed"] = pubmed_refs
    else:
        evidence["warnings"].append("No PubMed references available.")
        
    # Calculate score
    score, strength = calculate_evidence_score(evidence)
    evidence["evidence_score"] = score
    evidence["evidence_strength"] = strength
    
    return evidence
