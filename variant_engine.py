"""
variant_engine.py — Core Bioinformatics Analysis Engine
=========================================================
Implements the analytical logic for SNP / CNV classification,
functional impact prediction, and report generation.

Biological foundations:
  - The Central Dogma: DNA → RNA → Protein
  - A SNP changes one nucleotide; its consequence depends on codon context.
  - Synonymous: codon still encodes the same amino acid (silent)
  - Missense  : codon encodes a different amino acid (possibly harmful)
  - Nonsense  : codon becomes a premature stop (usually loss-of-function)
  - CNVs (Copy Number Variants): regions duplicated or deleted;
    dosage of the encoded genes is altered.
"""

import hashlib
from db_handler import (
    get_dbsnp_record,
    get_clinvar_record,
    get_gene_info,
    map_position_to_gene,
    get_chromosome_info,
)

# ── Codon / amino acid look-up tables ─────────────────────────────────────────
# Standard genetic code: each 3-base codon → 1-letter amino acid
# Used to deterministically "predict" synonymous vs missense vs nonsense
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# ── Clinical significance badge colours ───────────────────────────────────────
SIGNIFICANCE_COLOURS = {
    "Pathogenic":            "#e74c3c",
    "Likely_pathogenic":     "#e67e22",
    "Uncertain_significance":"#f39c12",
    "Likely_benign":         "#2ecc71",
    "Benign":                "#27ae60",
    "Risk_factor":           "#9b59b6",
    "drug_response":         "#3498db",
}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. VARIANT CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════════════════

def classify_variant(chromosome: str, position: int,
                     ref: str, alt: str,
                     cnv_type: str = None) -> dict:
    """
    Determine whether an input is a SNP or a CNV and return basic metadata.

    Rules:
      - If cnv_type is provided → CNV
      - If len(ref) == 1 and len(alt) == 1 → SNP (single nucleotide polymorphism)
      - If len(ref) != len(alt) → Indel (treated as small structural variant)
      - Otherwise → complex variant

    Returns a dict with keys:
        type        : "SNP" | "CNV" | "Indel" | "Complex"
        subtype     : "Transition" | "Transversion" | "Duplication" | "Deletion" | ...
        length      : 1 for SNP, bp span for CNV
    """
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    if cnv_type:
        return {
            "type": "CNV",
            "subtype": cnv_type.capitalize(),   # "Duplication" or "Deletion"
            "length": abs(len(ref) - len(alt)) if ref and alt else None,
        }

    ref, alt = ref.upper(), alt.upper()

    if len(ref) == 1 and len(alt) == 1:
        pair = (ref, alt)
        subtype = "Transition" if pair in transitions else "Transversion"
        return {"type": "SNP", "subtype": subtype, "length": 1}

    if len(ref) != len(alt):
        subtype = "Insertion" if len(alt) > len(ref) else "Deletion"
        return {"type": "Indel", "subtype": subtype,
                "length": abs(len(alt) - len(ref))}

    return {"type": "Complex", "subtype": "MNP", "length": len(ref)}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. FUNCTIONAL CONSEQUENCE PREDICTION
# ═══════════════════════════════════════════════════════════════════════════════

def predict_functional_impact(ref: str, alt: str,
                               consequence_hint: str = None) -> dict:
    """
    Predict the functional consequence of an SNP at the protein level.

    If a consequence_hint is provided (from dbSNP annotation), we trust that
    directly and enrich it with biological description.
    Otherwise we use a deterministic codon-simulation approach:
      1. Build a synthetic reference codon using position hash (for demo).
      2. Substitute the alternate base.
      3. Compare translated amino acids.

    Returns dict with:
        consequence : "synonymous_variant" | "missense_variant" |
                      "stop_gained" | "frameshift_variant" | "intron_variant" etc.
        impact_level: "LOW" | "MODERATE" | "HIGH" | "MODIFIER"
        description : Plain-English explanation
        sift_pred   : Simulated SIFT prediction ("Tolerated" / "Damaging")
        polyphen_pred: Simulated PolyPhen-2 prediction
    """

    # If dbSNP already tells us the consequence, map it directly
    consequence_map = {
        "synonymous_variant": {
            "impact_level": "LOW",
            "description": (
                "The nucleotide substitution results in the same amino acid "
                "due to codon degeneracy (silent mutation). Protein structure "
                "and function are likely unaffected."
            ),
            "sift_pred": "Tolerated",
            "polyphen_pred": "Benign",
        },
        "missense_variant": {
            "impact_level": "MODERATE",
            "description": (
                "A single nucleotide change alters one amino acid in the "
                "protein sequence. Depending on the physicochemical properties "
                "of the substitution and the residue's structural role, this "
                "may disrupt protein folding, binding, or enzymatic activity."
            ),
            "sift_pred": "Damaging",
            "polyphen_pred": "Possibly_damaging",
        },
        "stop_gained": {
            "impact_level": "HIGH",
            "description": (
                "A premature stop codon (nonsense mutation) truncates the "
                "protein. The truncated product is typically unstable, degraded "
                "by nonsense-mediated mRNA decay (NMD), or produces a "
                "non-functional protein — usually a loss-of-function event."
            ),
            "sift_pred": "Damaging",
            "polyphen_pred": "Probably_damaging",
        },
        "frameshift_variant": {
            "impact_level": "HIGH",
            "description": (
                "An insertion or deletion shifts the reading frame, altering "
                "all downstream codons. This almost always results in a "
                "completely aberrant protein sequence and/or premature termination."
            ),
            "sift_pred": "Damaging",
            "polyphen_pred": "Probably_damaging",
        },
        "intron_variant": {
            "impact_level": "MODIFIER",
            "description": (
                "The variant lies within an intronic region. While most intronic "
                "variants are non-coding, some affect splicing branch points, "
                "enhancers, or regulatory elements, potentially altering mRNA "
                "splicing or gene expression levels."
            ),
            "sift_pred": "Tolerated",
            "polyphen_pred": "Benign",
        },
        "splice_region_variant": {
            "impact_level": "LOW",
            "description": (
                "Located near an exon–intron boundary; may disrupt canonical "
                "splice signals (GT/AG rule), leading to exon skipping, "
                "intron retention, or activation of cryptic splice sites."
            ),
            "sift_pred": "Tolerated",
            "polyphen_pred": "Possibly_damaging",
        },
    }

    if consequence_hint and consequence_hint in consequence_map:
        result = dict(consequence_map[consequence_hint])
        result["consequence"] = consequence_hint
        return result

    # ── Fallback: simulate codon analysis ─────────────────────────────────────
    ref = ref.upper()
    alt = alt.upper()

    # Synthetic reference codon — first base is the variant position
    # We deterministically pick the other two bases using a simple hash
    codon_bases = "ACGT"
    seed = int(hashlib.md5(f"{ref}{alt}".encode()).hexdigest(), 16)
    b2 = codon_bases[seed % 4]
    b3 = codon_bases[(seed // 4) % 4]

    ref_codon = ref + b2 + b3
    alt_codon = alt + b2 + b3

    ref_aa = _CODON_TABLE.get(ref_codon, "?")
    alt_aa = _CODON_TABLE.get(alt_codon, "?")

    if alt_aa == "*":
        consequence = "stop_gained"
    elif ref_aa == alt_aa:
        consequence = "synonymous_variant"
    else:
        consequence = "missense_variant"

    result = dict(consequence_map[consequence])
    result["consequence"] = consequence
    result["ref_codon"] = ref_codon
    result["alt_codon"] = alt_codon
    result["ref_aa"] = ref_aa
    result["alt_aa"] = alt_aa
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# 3. SNP ANNOTATION
# ═══════════════════════════════════════════════════════════════════════════════

def annotate_snp(chromosome: str, position: int,
                 ref: str, alt: str) -> dict:
    """
    Full annotation pipeline for a user-submitted SNP.

    Steps:
      1. Classify the variant (SNP / Indel / Complex)
      2. Map position to a gene (Ensembl-like lookup)
      3. Predict functional consequence
      4. Look up gene details
      5. Assemble the annotation result

    Returns a comprehensive annotation dict suitable for rendering.
    """
    # Step 1 – Classify
    classification = classify_variant(chromosome, position, ref, alt)

    # Step 2 – Map to gene
    gene_record = map_position_to_gene(chromosome, position)

    # Step 3 – Functional impact
    impact = predict_functional_impact(ref, alt)

    # Step 4 – Gene details (if mapped)
    gene_info = get_gene_info(gene_record["gene"]) if gene_record else None

    # Step 5 – Chromosome metadata
    chrom_info = get_chromosome_info(chromosome)

    return {
        "input": {
            "chromosome": chromosome,
            "position": position,
            "ref": ref.upper(),
            "alt": alt.upper(),
        },
        "classification": classification,
        "gene": gene_record,
        "gene_info": gene_info,
        "impact": impact,
        "chromosome_info": chrom_info,
        "source": "local_coordinates",
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 4. CNV ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_cnv(chromosome: str, start: int, end: int,
                cnv_type: str, copy_number: int = None) -> dict:
    """
    Analyse a Copy Number Variant (CNV) — a deletion or duplication spanning
    a genomic region.

    Biological background:
      - Normal diploid humans carry 2 copies of each autosomal gene.
      - Duplications (copy number > 2) can over-express gene products.
      - Deletions (copy number < 2) can reduce or eliminate gene expression.
      - Haploinsufficiency: one copy is insufficient for normal function.
      - Triplosensitivity: one extra copy disrupts normal dosage balance.

    Parameters:
        chromosome : e.g. "17"
        start      : genomic start coordinate
        end        : genomic end coordinate
        cnv_type   : "Duplication" or "Deletion"
        copy_number: integer (e.g., 3 for dup, 1 for het del, 0 for hom del)

    Returns a detailed CNV annotation dict.
    """
    region_size = end - start
    gene_record = map_position_to_gene(chromosome, (start + end) // 2)

    # Dosage effect prediction
    if cnv_type.lower() == "duplication":
        if copy_number and copy_number >= 4:
            dosage_effect = "Severe triplosensitivity — likely pathogenic overexpression"
            dosage_class  = "HIGH"
        else:
            dosage_effect = "Gene dosage increase — possible overexpression phenotype"
            dosage_class  = "MODERATE"
    else:  # deletion
        if copy_number == 0:
            dosage_effect = "Complete gene loss (homozygous deletion) — high pathogenicity risk"
            dosage_class  = "HIGH"
        elif copy_number == 1:
            dosage_effect = "Heterozygous deletion — haploinsufficiency risk if gene is dosage-sensitive"
            dosage_class  = "MODERATE"
        else:
            dosage_effect = "Partial deletion — functional impact depends on exon involvement"
            dosage_class  = "LOW"

    # Region size classification (used in real CNV reporting)
    if region_size < 1000:
        size_class = "Micro-CNV (<1 kb)"
    elif region_size < 50000:
        size_class = "Small CNV (1–50 kb)"
    elif region_size < 500000:
        size_class = "Medium CNV (50–500 kb)"
    else:
        size_class = "Large CNV (>500 kb)"

    return {
        "input": {
            "chromosome": chromosome,
            "start": start,
            "end": end,
            "cnv_type": cnv_type,
            "copy_number": copy_number,
        },
        "region_size_bp": region_size,
        "size_class": size_class,
        "gene": gene_record,
        "dosage_effect": dosage_effect,
        "dosage_class": dosage_class,
        "classification": {
            "type": "CNV",
            "subtype": cnv_type.capitalize(),
        },
        "clinical_note": (
            f"A {region_size:,} bp {cnv_type.lower()} on chromosome {chromosome} "
            f"({'encompassing ' + gene_record['gene'] if gene_record else 'intergenic region'}). "
            f"Copy number inferred: {copy_number if copy_number is not None else 'not specified'}."
        ),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 5. BATCH rsID ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def batch_analyze_rsids(rsid_list: list[str]) -> list[dict]:
    """
    Process a list of rsIDs through the annotation pipeline.

    For each rsID:
      1. Look up dbSNP record → get basic variant info
      2. Look up ClinVar record → get clinical significance
      3. Predict functional impact (using consequence from dbSNP)
      4. Retrieve gene info

    Returns a list of result dicts (one per rsID), suitable for table display.
    """
    results = []

    for rsid in rsid_list:
        rsid = rsid.strip()
        if not rsid:
            continue

        snp_record    = get_dbsnp_record(rsid)
        clinvar_record = get_clinvar_record(rsid)

        if snp_record:
            impact = predict_functional_impact(
                snp_record.get("ref", "A"),
                snp_record.get("alt", "G"),
                snp_record.get("consequence"),
            )
            gene_info = get_gene_info(snp_record.get("gene", ""))

            results.append({
                "rsid": rsid,
                "found": True,
                "chromosome": snp_record.get("chromosome"),
                "position": snp_record.get("position"),
                "ref": snp_record.get("ref"),
                "alt": snp_record.get("alt"),
                "gene": snp_record.get("gene"),
                "consequence": snp_record.get("consequence", "unknown"),
                "hgvs": snp_record.get("hgvs"),
                "protein_change": snp_record.get("protein_change"),
                "allele_frequency": snp_record.get("allele_frequency"),
                "impact_level": impact.get("impact_level", "UNKNOWN"),
                "clinical_significance": (
                    clinvar_record.get("clinical_significance") if clinvar_record else "Not in ClinVar"
                ),
                "condition": (
                    clinvar_record.get("condition") if clinvar_record else "—"
                ),
                "gene_description": (
                    gene_info.get("description") if gene_info else "—"
                ),
                "pathway": (
                    gene_info.get("pathway") if gene_info else "—"
                ),
                "review_status": (
                    clinvar_record.get("review_status") if clinvar_record else "—"
                ),
                "functional_consequence_note": impact.get("description", ""),
            })
        else:
            results.append({
                "rsid": rsid,
                "found": False,
                "chromosome": None,
                "position": None,
                "ref": None,
                "alt": None,
                "gene": "—",
                "consequence": "Not found",
                "hgvs": "—",
                "protein_change": "—",
                "allele_frequency": None,
                "impact_level": "UNKNOWN",
                "clinical_significance": "Not in database",
                "condition": "—",
                "gene_description": "—",
                "pathway": "—",
                "review_status": "—",
                "functional_consequence_note": (
                    f"rsID '{rsid}' was not found in the local mock dataset. "
                    "In a production system, this query would be forwarded to "
                    "the NCBI E-utilities API."
                ),
            })

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# 6. FULL REPORT GENERATION (rsID-based)
# ═══════════════════════════════════════════════════════════════════════════════

def generate_interpretation_report(rsid: str) -> dict:
    """
    Generate a comprehensive variant interpretation report for a single rsID.

    Combines dbSNP + ClinVar + gene coordinates into a unified report dict
    that powers the detailed results page.
    """
    snp_record     = get_dbsnp_record(rsid)
    clinvar_record = get_clinvar_record(rsid)

    if not snp_record:
        return {"found": False, "rsid": rsid}

    impact    = predict_functional_impact(
        snp_record.get("ref", "A"),
        snp_record.get("alt", "G"),
        snp_record.get("consequence"),
    )
    gene_info = get_gene_info(snp_record.get("gene", ""))
    classification = classify_variant(
        snp_record.get("chromosome", "?"),
        snp_record.get("position", 0),
        snp_record.get("ref", "A"),
        snp_record.get("alt", "G"),
    )

    significance = (
        clinvar_record.get("clinical_significance") if clinvar_record
        else "Not assessed"
    )
    sig_colour = SIGNIFICANCE_COLOURS.get(significance, "#95a5a6")

    # ── Build a scientific interpretation paragraph ────────────────────────────
    gene_name = snp_record.get("gene", "unknown gene")
    consequence = snp_record.get("consequence", "unknown consequence")
    protein_change = snp_record.get("protein_change") or "N/A"
    hgvs = snp_record.get("hgvs") or "N/A"
    freq = snp_record.get("allele_frequency")
    freq_str = f"{freq:.4f} ({freq*100:.2f}%)" if freq is not None else "unknown"

    interpretation = (
        f"Variant {rsid} ({hgvs}) introduces a {classification['subtype'].lower()} "
        f"in {gene_name}, resulting in a {consequence.replace('_', ' ')} "
        f"({protein_change}). "
        f"The global alternate allele frequency is {freq_str}. "
    )
    if clinvar_record:
        interpretation += (
            f"ClinVar classifies this variant as '{significance}' in the context of "
            f"{clinvar_record.get('condition', 'unspecified condition')}. "
            f"The functional consequence is reported as: "
            f"{clinvar_record.get('functional_consequence', 'not described')}."
        )
    if gene_info:
        interpretation += (
            f" {gene_name} encodes {gene_info.get('description', '')}. "
            f"It participates in the {gene_info.get('pathway', 'unknown')} pathway."
        )

    return {
        "found": True,
        "rsid": rsid,
        "snp": snp_record,
        "clinvar": clinvar_record,
        "gene_info": gene_info,
        "classification": classification,
        "impact": impact,
        "clinical_significance": significance,
        "sig_colour": sig_colour,
        "interpretation": interpretation,
    }
