"""
db_handler.py — Mock Database Handler
======================================
Simulates queries to dbSNP, ClinVar, and Ensembl-style gene coordinate data.
All data is loaded from local JSON files. The architecture is designed so that
real API calls (e.g., NCBI E-utilities, Ensembl REST) can replace these functions
without changing the calling interface in variant_engine.py.

Bioinformatics context:
  - dbSNP: NCBI's database of Short Genetic Variations (SNPs, indels)
  - ClinVar: NCBI's archive of variant-disease relationships
  - Ensembl: Genome annotation and gene coordinate database
"""

import json
import os

# ── Path configuration ─────────────────────────────────────────────────────────
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
DBSNP_PATH   = os.path.join(DATA_DIR, "dbsnp_mock.json")
CLINVAR_PATH = os.path.join(DATA_DIR, "clinvar_mock.json")
GENES_PATH   = os.path.join(DATA_DIR, "gene_coordinates.json")


def _load_json(path: str) -> dict:
    """Utility: load a JSON file and return its contents as a dict."""
    with open(path, "r") as fh:
        return json.load(fh)


# ── Cached data loads (loaded once at startup) ─────────────────────────────────
_DBSNP_DATA   = _load_json(DBSNP_PATH)
_CLINVAR_DATA = _load_json(CLINVAR_PATH)
_GENE_DATA    = _load_json(GENES_PATH)


# ═══════════════════════════════════════════════════════════════════════════════
# dbSNP queries
# ═══════════════════════════════════════════════════════════════════════════════

def get_dbsnp_record(rsid: str) -> dict | None:
    """
    Retrieve a variant record from the mock dbSNP dataset by rsID.

    Real-world equivalent:
        GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
            ?db=snp&id=<rsid>&rettype=json

    Returns:
        dict with keys: rsid, chromosome, position, ref, alt, gene,
                        consequence, hgvs, protein_change, allele_frequency
        None if rsID not found.
    """
    # Normalise: strip "rs" prefix if user typed "rs1234", keep lowercase
    key = rsid.lower().strip()
    if not key.startswith("rs"):
        key = "rs" + key

    return _DBSNP_DATA.get(key)


def list_all_rsids() -> list[str]:
    """Return all rsIDs available in the mock dbSNP dataset."""
    return list(_DBSNP_DATA.keys())


# ═══════════════════════════════════════════════════════════════════════════════
# ClinVar queries
# ═══════════════════════════════════════════════════════════════════════════════

def get_clinvar_record(rsid: str) -> dict | None:
    """
    Retrieve a clinical significance record by rsID from mock ClinVar data.

    Real-world equivalent:
        GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi
            ?db=clinvar&term=<rsid>[VARID]

    Returns:
        dict with keys: clinical_significance, review_status, condition,
                        condition_id, last_evaluated, accession
        None if rsID not found.
    """
    key = rsid.lower().strip()
    if not key.startswith("rs"):
        key = "rs" + key

    return _CLINVAR_DATA.get(key)


# ═══════════════════════════════════════════════════════════════════════════════
# Gene coordinate / Ensembl-style queries
# ═══════════════════════════════════════════════════════════════════════════════

def get_gene_info(gene_symbol: str) -> dict | None:
    """
    Retrieve genomic coordinates and functional annotation for a gene symbol.

    Real-world equivalent:
        GET https://rest.ensembl.org/lookup/symbol/homo_sapiens/<gene>
            ?content-type=application/json

    Returns:
        dict with keys: gene, chromosome, start, end, strand,
                        biotype, description, pathway, omim
        None if gene not in local dataset.
    """
    symbol = gene_symbol.upper().strip()
    for entry in _GENE_DATA["gene_map"]:
        if entry["gene"].upper() == symbol:
            return entry
    return None


def map_position_to_gene(chromosome: str, position: int) -> dict | None:
    """
    Given a chromosomal position, find the overlapping gene in our dataset.
    Checks if position falls within [gene.start, gene.end] on the same chromosome.

    This is a simplified linear scan. Real tools use interval trees (e.g., PyRanges)
    for O(log n) lookup over all ~20,000 human genes.

    Returns:
        Gene info dict if position overlaps a known gene, else None.
    """
    chrom = chromosome.upper().replace("CHR", "")
    for entry in _GENE_DATA["gene_map"]:
        if entry["chromosome"].upper() == chrom:
            if entry["start"] <= position <= entry["end"]:
                return entry
    return None


def get_chromosome_info(chromosome: str) -> dict | None:
    """Return metadata (length, gene count) for a chromosome."""
    chrom = chromosome.upper().replace("CHR", "")
    return _GENE_DATA["chromosomes"].get(chrom)


# ═══════════════════════════════════════════════════════════════════════════════
# Bulk / batch helpers
# ═══════════════════════════════════════════════════════════════════════════════

def get_all_known_rsids() -> list[str]:
    """Return sorted list of all rsIDs in mock dataset (for demo / autocomplete)."""
    return sorted(_DBSNP_DATA.keys())
