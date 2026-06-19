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
import re
import requests
import xml.etree.ElementTree as ET
import urllib.parse
import time
import logging
from functools import lru_cache

GENE_FALLBACKS = {
    "rs6025": "F5",
    "rs7903146": "TCF7L2"
}

logger = logging.getLogger(__name__)

# ── Path configuration ─────────────────────────────────────────────────────────
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
DBSNP_PATH   = os.path.join(DATA_DIR, "dbsnp_mock.json")
CLINVAR_PATH = os.path.join(DATA_DIR, "clinvar_mock.json")
GENES_PATH   = os.path.join(DATA_DIR, "gene_coordinates.json")


def _load_json(path: str) -> dict:
    """Utility: load a JSON file and return its contents as a dict."""
    with open(path, "r") as fh:
        return json.load(fh)

def fetch_with_backoff(url: str, source: str, rsid: str, max_retries: int = 3, timeout: int = 10, response_type: str = "json"):
    """Fetch from URL with exponential backoff for rate limits."""
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=timeout)
            if resp.status_code == 429:
                logger.warning(f"{source} rate limit encountered for {rsid}. Retry {attempt + 1}/{max_retries}.")
                time.sleep((2 ** attempt) * 0.5)
                continue
            resp.raise_for_status()
            if response_type == "json":
                return resp.json()
            return resp.text
        except (requests.RequestException, json.JSONDecodeError) as e:
            if attempt == max_retries - 1:
                return None
            time.sleep((2 ** attempt) * 0.5)
    return None


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

    # 1. Check local mock dataset first
    record = _DBSNP_DATA.get(key)
    if record:
        # Merge live ClinVar data into mock record
        record["clinvar"] = _fetch_live_clinvar(key)
        return record
    
    # 2. Fallback to live NCBI Variation Services (Phase 1)
    return _fetch_live_dbsnp(key)

@lru_cache(maxsize=128)
def _fetch_live_dbsnp(rsid: str) -> dict | None:
    """
    Fetch variant information from NCBI E-utilities (esummary) for a given rsID.
    Parses the JSON payload to match the standard db_handler schema.
    """
    rsid_num = rsid.lower().replace("rs", "")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=snp&id={rsid_num}&retmode=json"
    
    data_json = fetch_with_backoff(url, "NCBI dbSNP", rsid)
    if not data_json:
        return None
        
    data = data_json.get("result", {}).get(rsid_num)
    if not data or "error" in data:
        return None

    # Parse basic fields
    chromosome = data.get("chr", "")
    chrpos = data.get("chrpos", "")
    position = int(chrpos.split(":")[1]) if ":" in chrpos else None

    # Parse alleles from SPDI
    ref, alt = "", ""
    spdi = data.get("spdi", "")
    if spdi:
        parts = spdi.split(",")[0].split(":")
        if len(parts) >= 4:
            ref = parts[2]
            alt = parts[3]

    # Parse docsum for HGVS, Protein change, and regex fallbacks
    docsum = data.get("docsum", "")

    # Multi-tiered Gene Extraction
    gene = ""
    # Tier 1: Extract from top-level "genes" array
    genes = data.get("genes", [])
    if genes:
        gene = genes[0].get("name", "")
    
    # Tier 2: Extract from primary_snapshot_data allele_annotations
    if not gene:
        try:
            snapshot = data.get("primary_snapshot_data", {})
            annotations = snapshot.get("allele_annotations", [])
            for ann in annotations:
                assembly_ann = ann.get("assembly_annotation", [])
                for aa in assembly_ann:
                    for g in aa.get("genes", []):
                        if g.get("symbol"):
                            gene = g.get("symbol")
                            break
                    if gene: break
                if gene: break
        except Exception:
            pass

    # Tier 3: Regex match inside docsum
    if not gene and docsum:
        match_gene = re.search(r"GENE=([^|]+)", docsum)
        if match_gene:
            gene_info_str = match_gene.group(1)
            gene = gene_info_str.split(":")[0] if ":" in gene_info_str else gene_info_str

    # Tier 4: Scientific Fallback Layer
    if not gene and f"rs{rsid_num}" in GENE_FALLBACKS:
        gene = GENE_FALLBACKS[f"rs{rsid_num}"]

    # Parse consequence
    consequence = data.get("fxn_class", "")

    hgvs = ""
    protein_change = ""
    if docsum:
        # Try to find HGVS NM_ or NC_
        match_hgvs = re.search(r"HGVS=([^|]+)", docsum)
        if match_hgvs:
            hgvs_list = match_hgvs.group(1).split(",")
            # prefer NM_ (transcript) over genomic coords if multiple exist
            for h in hgvs_list:
                if h.startswith("NM_"):
                    hgvs = h
                    break
            if not hgvs and hgvs_list:
                hgvs = hgvs_list[0]

        # Try to find protein change NP_
        match_p = re.search(r"(NP_[^|]+)", docsum)
        if match_p:
            protein_change = match_p.group(1).split(",")[0]
            if ":" in protein_change:
                protein_change = protein_change.split(":")[1]

    # Parse allele frequency
    global_mafs = data.get("global_mafs", [])
    allele_frequency = None
    if global_mafs:
        # Prefer major population cohorts
        freq_str = ""
        for study_name in ["GnomAD_genomes", "1000Genomes", "ExAC"]:
            for maf in global_mafs:
                if maf.get("study") == study_name:
                    freq_str = maf.get("freq", "")
                    break
            if freq_str:
                break
        
        # Fallback to the first available frequency
        if not freq_str and global_mafs:
            freq_str = global_mafs[0].get("freq", "")
        
        if freq_str and "=" in freq_str:
            freq_val = freq_str.split("=")[1].split("/")[0]
            try:
                allele_frequency = float(freq_val)
            except ValueError:
                pass

    pop_data = _fetch_live_population_frequencies(rsid)
    population_frequencies = pop_data["frequencies"]
    frequency_source = pop_data["source"]

    # Only fallback to allele_frequency if global frequency is missing AND Ensembl actually returned data
    if frequency_source != "Unavailable" and population_frequencies.get("global") is None and allele_frequency is not None:
        population_frequencies["global"] = allele_frequency

    gwas_associations = _fetch_live_gwas(rsid)
    clinvar_data = _fetch_live_clinvar(rsid)

    return {
        "rsid": f"rs{rsid_num}",
        "chromosome": chromosome,
        "position": position,
        "ref": ref,
        "alt": alt,
        "gene": gene,
        "consequence": consequence,
        "hgvs": hgvs,
        "protein_change": protein_change,
        "allele_frequency": allele_frequency,
        "population_frequencies": population_frequencies,
        "frequency_source": frequency_source,
        "gwas_associations": gwas_associations,
        "clinvar": clinvar_data
    }

@lru_cache(maxsize=128)
def _fetch_live_population_frequencies(rsid: str) -> dict:
    """
    Fetch population frequencies from Ensembl Variation REST API.
    Maps gnomAD and 1000Genomes data to canonical ethnic categories.
    """
    rsid_num = rsid.lower().replace("rs", "")
    url = f"https://rest.ensembl.org/variation/human/rs{rsid_num}?content-type=application/json;pops=1"
    
    pops_map = {
        "south_asian": ["gnomADe:sas", "gnomADg:sas", "1000GENOMES:phase_3:SAS"],
        "european": ["gnomADe:nfe", "gnomADg:nfe", "1000GENOMES:phase_3:EUR"],
        "east_asian": ["gnomADe:eas", "gnomADg:eas", "1000GENOMES:phase_3:EAS"],
        "african": ["gnomADe:afr", "gnomADg:afr", "1000GENOMES:phase_3:AFR"],
        "american": ["gnomADe:amr", "gnomADg:amr", "1000GENOMES:phase_3:AMR"]
    }
    
    result = {
        "frequencies": {
            "global": None,
            "south_asian": None,
            "european": None,
            "east_asian": None,
            "african": None,
            "american": None
        },
        "source": "Ensembl REST API (gnomAD / 1000Genomes)"
    }
    
    unavailable_result = {
        "frequencies": {
            "available": False,
            "error": "Ensembl lookup failed"
        },
        "source": "Unavailable"
    }
    
    try:
        resp = requests.get(url, timeout=10)
        if not resp.ok: 
            return unavailable_result
        data = resp.json()
        
        # Try to get global from MAF
        maf = data.get("MAF")
        if maf is not None:
            try:
                result["frequencies"]["global"] = float(maf)
            except ValueError:
                pass
            
        populations = data.get("populations", [])
        if not populations:
            return result
            
        # Map frequencies
        for target_pop, source_pops in pops_map.items():
            for sp in source_pops:
                matched = [p for p in populations if p.get("population") == sp and p.get("frequency") is not None]
                if matched:
                    # Take minor allele frequency
                    minor_freqs = [float(p["frequency"]) for p in matched if float(p.get("frequency", 0)) <= 0.5]
                    if minor_freqs:
                        result["frequencies"][target_pop] = max(minor_freqs)
                    break
                    
    except (requests.RequestException, json.JSONDecodeError):
        return unavailable_result
        
    return result


@lru_cache(maxsize=128)
def _fetch_live_gwas(rsid: str) -> list[dict]:
    """
    Fetch quantitative association metrics from the GWAS Catalog REST API.
    Retrieves all associations, sorts them by p-value (strongest first),
    and fetches study details for the strongest association.
    """
    rsid_lower = rsid.lower()
    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid_lower}/associations"
    
    try:
        resp = requests.get(url, timeout=10)
        if not resp.ok:
            return []
            
        data = resp.json()
        assocs = data.get("_embedded", {}).get("associations", [])
        if not assocs:
            return []
            
        parsed_assocs = []
        for assoc in assocs:
            # We must have a pvalue to evaluate strength
            if assoc.get("pvalue") is None:
                continue
                
            item = {
                "p_value": assoc.get("pvalue"),
                "beta": assoc.get("betaNum"),
                "odds_ratio": assoc.get("orPerCopyNum"),
                "beta_unit": assoc.get("betaUnit"),
                "beta_direction": assoc.get("betaDirection"),
                "trait": None,
                "study": None,
                "pmid": None,
                "_study_url": None
            }
            
            # Extract study URL
            links = assoc.get("_links", {})
            study_link = links.get("study", {}).get("href")
            if study_link:
                # Strip templated query params if present
                if "{?projection}" in study_link:
                    study_link = study_link.replace("{?projection}", "")
                item["_study_url"] = study_link
                
            parsed_assocs.append(item)
            
        # Sort by p-value (lowest first)
        parsed_assocs.sort(key=lambda x: x["p_value"])
        
        # Only fetch study details for the top association to save time
        if parsed_assocs and parsed_assocs[0]["_study_url"]:
            top = parsed_assocs[0]
            try:
                study_resp = requests.get(top["_study_url"], timeout=10)
                if study_resp.ok:
                    study_data = study_resp.json()
                    
                    # Trait is typically in diseaseTrait.trait
                    trait = study_data.get("diseaseTrait", {}).get("trait")
                    pub_info = study_data.get("publicationInfo", {})
                    
                    top["trait"] = trait
                    top["study"] = pub_info.get("title")
                    top["pmid"] = pub_info.get("pubmedId")
            except (requests.RequestException, json.JSONDecodeError):
                pass
                
        # Cleanup internal fields
        for item in parsed_assocs:
            item.pop("_study_url", None)
            
        return parsed_assocs
        
    except (requests.RequestException, json.JSONDecodeError):
        return []


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


# ═══════════════════════════════════════════════════════════════════════════════
# Mock Evidence APIs
# ═══════════════════════════════════════════════════════════════════════════════

def get_pubmed_mock(rsid: str, gene: str) -> list[dict]:
    """
    Simulate fetching PubMed references for a variant or gene.
    Returns realistic mocked publications instead of 404 links.
    """
    key = rsid.lower().strip() if rsid else ""
    if not key.startswith("rs") and key:
        key = "rs" + key

    # Return predefined mock publications to prevent 404s
    pubmed_data = {
        "rs429358": [
            {
                "pmid": "33340485",
                "id": "PMID: 33340485",
                "title": "APOE and Alzheimer's disease: advances in genetics, pathophysiology, and therapeutic approaches.",
                "year": 2021,
                "url": "https://pubmed.ncbi.nlm.nih.gov/33340485/"
            },
            {
                "pmid": "31367008",
                "id": "PMID: 31367008",
                "title": "Apolipoprotein E and Alzheimer disease: pathobiology and targeting strategies.",
                "year": 2019,
                "url": "https://pubmed.ncbi.nlm.nih.gov/31367008/"
            },
            {
                "pmid": "30776012",
                "id": "PMID: 30776012",
                "title": "Is Alzheimer's Disease Risk Modifiable?",
                "year": 2019,
                "url": "https://pubmed.ncbi.nlm.nih.gov/30776012/"
            }
        ],
        "rs334": [
            {
                "pmid": "29542687",
                "id": "PMID: 29542687",
                "title": "Sickle cell disease.",
                "year": 2018,
                "url": "https://pubmed.ncbi.nlm.nih.gov/29542687/"
            },
            {
                "pmid": "28159390",
                "id": "PMID: 28159390",
                "title": "Sickle cell disease.",
                "year": 2017,
                "url": "https://pubmed.ncbi.nlm.nih.gov/28159390/"
            }
        ],
        "rs1800562": [
            {
                "pmid": "20301613",
                "id": "PMID: 20301613",
                "title": "HFE-Related Hemochromatosis.",
                "year": 2000,
                "url": "https://pubmed.ncbi.nlm.nih.gov/20301613/"
            },
            {
                "pmid": "30798813",
                "id": "PMID: 30798813",
                "title": "Hemochromatosis: Hereditary hemochromatosis and HFE gene.",
                "year": 2019,
                "url": "https://pubmed.ncbi.nlm.nih.gov/30798813/"
            }
        ]
    }
    
    # Return matched references, or empty list if no exact mock exists
    # This prevents the creation of broken links based on hashes
    return pubmed_data.get(key, [])

def get_genecards_mock(gene: str) -> dict | None:
    """
    Simulate fetching GeneCards summary for a gene.
    """
    if not gene or gene == "—":
        return None
        
    symbol = gene.upper().strip()
    
    # Enrich the gene coordinates data we already have
    gene_info = get_gene_info(symbol)
    if not gene_info:
        return None
        
    # Build a GeneCards-style response
    return {
        "id": f"{symbol} GeneCards Entry",
        "url": f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={symbol}",
        "gene_summary": gene_info.get("description", ""),
        "biological_function": gene_info.get("pathway", ""),
        "associated_diseases": [] # Normally extracted from GeneCards/MalaCards
    }

@lru_cache(maxsize=128)
def _fetch_live_clinvar(rsid: str) -> dict | None:
    """
    Fetches live ClinVar evidence via NCBI E-Utilities using ESearch, ESummary, and EFetch (XML).
    """
    try:
        # 1. ESearch
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={rsid}[rs]&retmode=json"
        search_json = fetch_with_backoff(search_url, "NCBI ClinVar (ESearch)", rsid)
        if not search_json:
            return None
            
        id_list = search_json.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return None
            
        # 2. ESummary
        summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={','.join(id_list)}&retmode=json"
        summary_json = fetch_with_backoff(summary_url, "NCBI ClinVar (ESummary)", rsid)
        if not summary_json:
            return None
            
        result = summary_json.get("result", {})
        
        target_uid = None
        target_data = None
        
        for uid in id_list:
            data = result.get(uid, {})
            if data.get("obj_type") == "single nucleotide variant":
                target_uid = uid
                target_data = data
                break
                
        # Fallback to haplotype or first result if no single nucleotide variant is explicitly found
        if not target_uid and id_list:
            target_uid = id_list[0]
            target_data = result.get(target_uid, {})
            
        if not target_uid or not target_data:
            return None
            
        gc = target_data.get("germline_classification", {})
        clinical_significance = gc.get("description", "Unknown")
        raw_review_status = gc.get("review_status", "no assertion")
        last_evaluated = gc.get("last_evaluated", "Unknown")
        
        traits = gc.get("trait_set", [])
        condition = traits[0].get("trait_name") if traits else "Not specified"
        
        accession = target_data.get("accession_version", "")
        
        # 3. EFetch (XML) for assertions
        efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid=1&id={target_uid}"
        fetch_text = fetch_with_backoff(efetch_url, "NCBI ClinVar (EFetch)", rsid, response_type="text")
        if not fetch_text:
            return None
            
        root = ET.fromstring(fetch_text)
        
        assertions = []
        for assertion_node in root.findall('.//ClinicalAssertion'):
            # Grab submitter name
            submitter = "Unknown"
            clinvar_accession = assertion_node.find('ClinVarAccession')
            if clinvar_accession is not None:
                submitter = clinvar_accession.get('SubmitterName', 'Unknown')
                
            # Grab classification
            classification_node = assertion_node.find('Classification')
            classification = "Unknown"
            if classification_node is not None:
                germ = classification_node.find('GermlineClassification')
                if germ is not None and germ.text:
                    classification = germ.text
                else:
                    clin_sig = classification_node.find('ClinicalSignificance')
                    if clin_sig is not None and clin_sig.text:
                        classification = clin_sig.text
            
            if classification != "Unknown":
                # Normalize case for aggregation
                norm_class = classification.capitalize()
                assertions.append({
                    "submitter": submitter,
                    "classification": norm_class
                })
                
        # 4. Aggregation and Scoring
        from collections import Counter
        class_counts = Counter(a["classification"] for a in assertions)
        
        total_assertions = sum(class_counts.values())
        submission_count = total_assertions
        
        consensus_score = None
        if total_assertions > 0:
            dominant_count = class_counts.most_common(1)[0][1]
            consensus_score = int(round((dominant_count / total_assertions) * 100))
            
        conflicting = "conflicting" in raw_review_status.lower()
        
        # Confidence Level Logic
        rl = raw_review_status.lower()
        if "practice guideline" in rl:
            confidence = "Very High"
        elif "expert panel" in rl:
            confidence = "High"
        elif "multiple submitters" in rl and not conflicting:
            confidence = "Moderate"
        else:
            confidence = "Low"
            
        return {
            "source": "live",
            "clinical_significance": clinical_significance,
            "review_status": raw_review_status,
            "raw_review_status": raw_review_status,
            "condition": condition,
            "accession": accession,
            "variation_id": target_uid,
            "submission_count": submission_count,
            "last_evaluated": last_evaluated,
            "confidence_level": confidence,
            "consensus_score": consensus_score,
            "conflicting_interpretations": conflicting,
            "conflict_breakdown": dict(class_counts) if conflicting or submission_count > 0 else {},
            "assertions": assertions
        }
    except Exception as e:
        print(f"Error fetching live ClinVar for {rsid}: {e}")
        return None

# ═══════════════════════════════════════════════════════════════════════════════
# PubMed queries
# ═══════════════════════════════════════════════════════════════════════════════

@lru_cache(maxsize=128)
def _fetch_live_pubmed(rsid: str, gene: str) -> dict:
    """
    Retrieve live literature evidence from PubMed combining rsID and gene symbol.
    Queries both to maximize literature retrieval, returning top 5 and total count.
    """
    fallback_response = {
        "available": False,
        "paper_count": 0,
        "papers": []
    }
    
    try:
        if not rsid and not gene:
            return fallback_response

        terms = []
        if rsid:
            terms.append(f"{rsid}[All Fields]")
        if gene:
            terms.append(f"{gene}[All Fields]")
            
        query_term = " OR ".join(terms)
        encoded_term = urllib.parse.quote(query_term)
        
        esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={encoded_term}&retmode=json&retmax=5"
        
        search_data = fetch_with_backoff(esearch_url, "NCBI PubMed ESearch", rsid, response_type="json")
        if not search_data:
            return fallback_response
            
        result = search_data.get("esearchresult", {})
        count_str = result.get("count", "0")
        try:
            total_count = int(count_str)
        except ValueError:
            total_count = 0
            
        idlist = result.get("idlist", [])
        
        if total_count == 0 or not idlist:
            return fallback_response
            
        # Fetch summaries for the top 5
        id_str = ",".join(idlist)
        esummary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={id_str}&retmode=json"
        
        summary_data = fetch_with_backoff(esummary_url, "NCBI PubMed ESummary", rsid, response_type="json")
        if not summary_data:
            return fallback_response
            
        papers = []
        summary_result = summary_data.get("result", {})
        
        for pmid in idlist:
            doc = summary_result.get(pmid, {})
            if not doc:
                continue
                
            title = doc.get("title", "Unknown Title")
            journal = doc.get("source", "Unknown Journal")
            
            pubdate = doc.get("pubdate", "")
            year = 0
            if pubdate:
                match = re.match(r'^(\d{4})', pubdate)
                if match:
                    year = int(match.group(1))
                    
            authors_list = doc.get("authors", [])
            authors = [a.get("name", "") for a in authors_list if a.get("name")]
            
            has_abstract = doc.get("hasabstract", 0) == 1
            
            papers.append({
                "pmid": pmid,
                "title": title,
                "journal": journal,
                "year": year,
                "authors": authors,
                "has_abstract": has_abstract,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            })
            
        return {
            "available": True,
            "paper_count": total_count,
            "papers": papers
        }
        
        return fallback_response
        
    except Exception as e:
        logger.error(f"Error fetching PubMed data for {rsid}/{gene}: {e}")
        return fallback_response

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 1: Gene Context
# ═══════════════════════════════════════════════════════════════════════════════

@lru_cache(maxsize=128)
def _fetch_live_gene_context(gene_symbol: str) -> dict:
    """
    Fetch gene context from NCBI Gene (Priority 1) or Ensembl (Priority 2).
    """
    fallback = {"available": False}
    if not gene_symbol:
        return fallback
        
    # Priority 1: NCBI Gene
    try:
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_symbol}[sym]+AND+human[orgn]&retmode=json"
        search_res = fetch_with_backoff(search_url, "NCBI Gene ESearch", gene_symbol)
        
        if search_res and "esearchresult" in search_res and search_res["esearchresult"].get("idlist"):
            gene_id = search_res["esearchresult"]["idlist"][0]
            summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json"
            summary_res = fetch_with_backoff(summary_url, "NCBI Gene ESummary", gene_symbol)
            
            if summary_res and "result" in summary_res and gene_id in summary_res["result"]:
                data = summary_res["result"][gene_id]
                return {
                    "available": True,
                    "symbol": data.get("nomenclaturesymbol", gene_symbol),
                    "full_name": data.get("nomenclaturename", ""),
                    "description": data.get("summary", ""),
                    "chromosome": str(data.get("chromosome", "")),
                    "location": data.get("maplocation", ""),
                    "organism": data.get("organism", {}).get("scientificname", "Homo sapiens"),
                    "source": "NCBI Gene"
                }
    except Exception as e:
        logger.error(f"NCBI Gene fetch failed for {gene_symbol}: {e}")

    # Priority 2: Ensembl
    try:
        ensembl_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1;content-type=application/json"
        resp = requests.get(ensembl_url, headers={"Content-Type": "application/json"}, timeout=10)
        if resp.ok:
            data = resp.json()
            return {
                "available": True,
                "symbol": data.get("display_name", gene_symbol),
                "full_name": data.get("description", "").split(" [")[0] if data.get("description") else "",
                "description": data.get("description", ""),
                "chromosome": str(data.get("seq_region_name", "")),
                "location": f"Chr{data.get('seq_region_name')}:{data.get('start')}-{data.get('end')}",
                "organism": "Homo sapiens",
                "source": "Ensembl"
            }
    except Exception as e:
        logger.error(f"Ensembl Gene fetch failed for {gene_symbol}: {e}")
        
    return fallback

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 2: Disease Association
# ═══════════════════════════════════════════════════════════════════════════════

@lru_cache(maxsize=1)
def _load_gene_diseases() -> dict:
    """Loads the external JSON database for gene diseases."""
    filepath = os.path.join(os.path.dirname(__file__), "data", "gene_diseases.json")
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to load gene_diseases.json: {e}")
        return {}

def _fetch_gene_diseases(gene_symbol: str) -> dict:
    """Retrieves disease associations for a gene from the knowledge base."""
    db = _load_gene_diseases()
    return db.get(gene_symbol, {"available": False})

# ═══════════════════════════════════════════════════════════════════════════════
# Phase 3: Biological Pathways
# ═══════════════════════════════════════════════════════════════════════════════

@lru_cache(maxsize=1)
def _load_gene_pathways() -> dict:
    """Loads the external JSON database for gene pathways."""
    filepath = os.path.join(os.path.dirname(__file__), "data", "gene_pathways.json")
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to load gene_pathways.json: {e}")
        return {}

def _fetch_gene_pathways(gene_symbol: str) -> dict:
    """Retrieves pathway associations for a gene from the knowledge base."""
    db = _load_gene_pathways()
    return db.get(gene_symbol, {"available": False})


_pgx_cache = None
def _load_pgx():
    global _pgx_cache
    if _pgx_cache is None:
        try:
            path = os.path.join(DATA_DIR, "pharmacogenomics.json")
            with open(path, "r", encoding="utf-8") as f:
                _pgx_cache = json.load(f)
        except Exception:
            _pgx_cache = {}
    return _pgx_cache

def _fetch_pharmacogenomics(gene_symbol: str) -> dict:
    if not gene_symbol:
        return {"available": False, "drug_gene_interactions": [], "applications": []}
        
    db = _load_pgx()
    record = db.get(gene_symbol)
    if record:
        return {
            "available": True,
            "drug_gene_interactions": record.get("drug_gene_interactions", []),
            "applications": record.get("applications", [])
        }
    return {"available": False, "drug_gene_interactions": [], "applications": []}


_panels_cache = None
def _load_panels():
    global _panels_cache
    if _panels_cache is None:
        try:
            path = os.path.join(DATA_DIR, "disease_panels.json")
            with open(path, "r", encoding="utf-8") as f:
                _panels_cache = json.load(f)
        except Exception:
            _panels_cache = {}
    return _panels_cache

def _fetch_disease_panel(disease_name: str) -> list:
    db = _load_panels()
    for k, v in db.items():
        if k.lower() == disease_name.lower():
            return v.get("genes", [])
    return []
