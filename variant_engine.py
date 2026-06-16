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
    _fetch_live_pubmed,
    _fetch_live_gene_context,
    _fetch_gene_diseases,
    _fetch_gene_pathways,
)
from evidence_manager import build_evidence_object

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

# ── Impact Explanation Reference ──────────────────────────────────────────────
# Reusable mapping that explains what each impact level means in clinical
# genomics context. Consumed by all analysis outputs (SNP, CNV, batch, report).
# This is explanatory metadata only — it does NOT modify interpretation logic.

IMPACT_EXPLANATIONS = {
    "HIGH": {
        "meaning": "Predicted to substantially alter or abolish gene/protein function.",
        "biological_consequences": [
            "Complete loss of protein product (nonsense-mediated decay)",
            "Severely truncated or non-functional protein",
            "Disruption of critical protein domains or active sites",
            "Loss of essential splice signals leading to aberrant mRNA",
        ],
        "example_variant_types": [
            "Stop gained (nonsense)",
            "Frameshift variant",
            "Splice donor/acceptor loss",
            "Start lost",
        ],
        "interpretation": (
            "Variants with HIGH impact have the strongest predicted effect on "
            "gene function. In clinical genomics, these are prioritised for "
            "further evaluation, especially when found in genes associated with "
            "Mendelian disease. A HIGH impact prediction alone does not confirm "
            "pathogenicity — clinical significance depends on the specific gene, "
            "inheritance pattern, and supporting evidence."
        ),
    },
    "MODERATE": {
        "meaning": "Predicted to alter protein sequence or function without complete disruption.",
        "biological_consequences": [
            "Single amino acid substitution (may affect folding or binding)",
            "In-frame insertion/deletion altering protein length",
            "Altered protein-protein interaction interfaces",
            "Changed enzymatic activity or substrate specificity",
        ],
        "example_variant_types": [
            "Missense variant",
            "In-frame deletion",
            "In-frame insertion",
            "Protein-altering variant",
        ],
        "interpretation": (
            "MODERATE impact variants change the protein but do not destroy it. "
            "The actual effect ranges from benign (tolerated substitution in a "
            "non-critical region) to damaging (substitution at an active site or "
            "conserved residue). Computational tools like SIFT and PolyPhen-2 "
            "help distinguish tolerated from deleterious changes, but functional "
            "studies provide the strongest evidence."
        ),
    },
    "LOW": {
        "meaning": "Predicted to have minimal or no effect on protein function.",
        "biological_consequences": [
            "No change to protein sequence (synonymous / silent)",
            "Variant in a region with low evolutionary conservation",
            "Possible minor effects on mRNA stability or splicing regulation",
            "Potential codon usage bias without phenotypic consequence",
        ],
        "example_variant_types": [
            "Synonymous variant (silent)",
            "Splice region variant (near but not at canonical site)",
            "5' UTR variant",
            "3' UTR variant",
        ],
        "interpretation": (
            "LOW impact variants are generally considered benign in clinical "
            "practice. While most are functionally neutral, rare exceptions "
            "exist — for example, synonymous variants that disrupt exonic "
            "splicing enhancers (ESEs) or alter mRNA secondary structure. "
            "These are typically not prioritised unless other evidence suggests "
            "a functional role."
        ),
    },
    "MODIFIER": {
        "meaning": "Variant in a non-coding or intergenic region with uncertain functional relevance.",
        "biological_consequences": [
            "No direct effect on protein-coding sequence",
            "Possible regulatory effects on gene expression",
            "Potential impact on enhancer, silencer, or insulator elements",
            "May affect non-coding RNA genes in the vicinity",
        ],
        "example_variant_types": [
            "Intron variant",
            "Intergenic variant",
            "Upstream/downstream gene variant",
            "Non-coding transcript variant",
        ],
        "interpretation": (
            "MODIFIER variants fall outside protein-coding regions and are "
            "the most common class of human genetic variation. The vast majority "
            "have no detectable phenotypic effect. However, a small subset may "
            "influence gene regulation — particularly variants in promoters, "
            "enhancers, or conserved non-coding elements identified by ENCODE "
            "or Roadmap Epigenomics data."
        ),
    },
}


def get_impact_explanation(impact_level: str) -> dict:
    """Return the explanation dict for a given impact level, with fallback."""
    return IMPACT_EXPLANATIONS.get(
        impact_level.upper() if impact_level else "MODIFIER",
        IMPACT_EXPLANATIONS["MODIFIER"],
    )


# ── Clinical Significance Explanation Reference ───────────────────────────────
# Reusable mapping that explains what each ClinVar/clinical classification means.
# Consumed by all analysis outputs (SNP, CNV, batch, report).
# This is explanatory metadata only — it does NOT modify classification logic.
#
# Future-compatible fields (currently null, populated when modules are added):
#   - acmg_criteria: ACMG/AMP rule codes (e.g. PS1, PM2, BP4)
#   - odds_ratio / beta_value / p_value: statistical effect-size metrics
#   - population_frequency_note: gnomAD/1000G context

SIGNIFICANCE_EXPLANATIONS = {
    "Pathogenic": {
        "definition": (
            "The variant has been determined to cause disease based on strong "
            "evidence from multiple independent sources."
        ),
        "interpretation": (
            "Pathogenic variants are considered causative for the associated "
            "condition. In clinical settings, a Pathogenic classification "
            "typically triggers diagnostic confirmation, cascade family testing, "
            "and disease-specific management. This is the highest-confidence "
            "classification for disease-causing variants."
        ),
        "evidence_basis": [
            "ClinVar expert panel or practice guideline review",
            "Functional studies demonstrating deleterious effect",
            "Segregation with disease in multiple affected families",
            "Absence or extremely low frequency in population databases (gnomAD)",
            "Published case-level evidence in peer-reviewed literature",
        ],
        "acmg_criteria": None,
        "odds_ratio": None,
        "beta_value": None,
        "p_value": None,
        "population_frequency_note": None,
    },
    "Likely_pathogenic": {
        "definition": (
            "The variant has strong but not definitive evidence of causing disease "
            "(typically >90% certainty)."
        ),
        "interpretation": (
            "Likely Pathogenic variants are treated similarly to Pathogenic in "
            "most clinical contexts. The distinction indicates that while evidence "
            "is compelling, additional data (e.g. functional studies, segregation "
            "in more families) could further strengthen or occasionally reclassify "
            "the finding. Periodic re-evaluation is recommended as new evidence "
            "emerges."
        ),
        "evidence_basis": [
            "ClinVar submitter consensus with multiple concordant assertions",
            "Computational predictions supporting deleteriousness (SIFT, PolyPhen-2)",
            "Variant located in a mutational hotspot or well-established functional domain",
            "Low allele frequency consistent with disease prevalence",
            "At least one published case report with phenotype correlation",
        ],
        "acmg_criteria": None,
        "odds_ratio": None,
        "beta_value": None,
        "p_value": None,
        "population_frequency_note": None,
    },
    "Uncertain_significance": {
        "definition": (
            "The evidence for or against pathogenicity is currently insufficient "
            "or conflicting — the clinical impact cannot be determined."
        ),
        "interpretation": (
            "Variants of Uncertain Significance (VUS) should NOT be used for "
            "clinical decision-making without additional evidence. VUS are common "
            "and reflect the current limits of knowledge. Many VUS are eventually "
            "reclassified as Benign or Pathogenic as evidence accumulates. "
            "Periodic re-evaluation through ClinVar or the testing laboratory "
            "is strongly recommended."
        ),
        "evidence_basis": [
            "Insufficient ClinVar submissions or conflicting interpretations",
            "Computational predictions are discordant or borderline",
            "Limited or no functional study data available",
            "Allele frequency data neither confirms nor excludes pathogenicity",
            "No published case reports establishing clear genotype-phenotype link",
        ],
        "acmg_criteria": None,
        "odds_ratio": None,
        "beta_value": None,
        "p_value": None,
        "population_frequency_note": None,
    },
    "Likely_benign": {
        "definition": (
            "The variant has strong evidence suggesting it does not cause disease "
            "(typically >90% certainty of benign status)."
        ),
        "interpretation": (
            "Likely Benign variants are generally not considered clinically "
            "actionable. They are common in population databases and lack "
            "evidence of functional disruption. However, the classification "
            "carries a small residual uncertainty, so the variant may warrant "
            "re-evaluation if new evidence emerges or if the clinical "
            "presentation is highly suggestive."
        ),
        "evidence_basis": [
            "Elevated allele frequency in population databases (gnomAD, 1000 Genomes)",
            "Computational predictions indicate tolerated/benign change",
            "No published disease associations in literature",
            "ClinVar submitters consistently classify as Benign/Likely Benign",
            "Functional studies show no detectable impact on protein function",
        ],
        "acmg_criteria": None,
        "odds_ratio": None,
        "beta_value": None,
        "p_value": None,
        "population_frequency_note": None,
    },
    "Benign": {
        "definition": (
            "The variant is established as not causing disease, supported by "
            "strong concordant evidence."
        ),
        "interpretation": (
            "Benign variants represent normal human genetic variation. They are "
            "not clinically actionable and should not influence medical management. "
            "Many Benign variants are common polymorphisms observed across "
            "diverse populations. This is the highest-confidence classification "
            "for non-pathogenic variants."
        ),
        "evidence_basis": [
            "High allele frequency in multiple population databases",
            "Multiple concordant Benign assertions in ClinVar",
            "No functional consequence predicted or demonstrated",
            "Observed in healthy individuals without the associated phenotype",
            "Well-characterized common polymorphism in literature",
        ],
        "acmg_criteria": None,
        "odds_ratio": None,
        "beta_value": None,
        "p_value": None,
        "population_frequency_note": None,
    },
    "Risk_factor": {
        "definition": (
            "The variant is statistically associated with increased susceptibility "
            "to a condition but is neither necessary nor sufficient to cause it."
        ),
        "interpretation": (
            "Risk Factor variants modify disease probability but do not "
            "deterministically cause disease. A carrier may never develop the "
            "associated condition, and many affected individuals do not carry the "
            "variant. Risk assessment requires integration with family history, "
            "environmental exposures, and polygenic risk scores. These variants "
            "are most useful in risk stratification, not diagnosis."
        ),
        "evidence_basis": [
            "Genome-wide association studies (GWAS) with significant p-values",
            "Replicated odds ratios / relative risk across multiple cohorts",
            "Population-level allele frequency differences between cases and controls",
            "Functional studies suggesting mechanistic plausibility",
            "ClinVar classification as Risk_factor with supporting evidence",
        ],
        "acmg_criteria": None,
        "odds_ratio": None,
        "beta_value": None,
        "p_value": None,
        "population_frequency_note": None,
    },
    "drug_response": {
        "definition": (
            "The variant is associated with altered pharmacological response to "
            "one or more medications — affecting drug efficacy, metabolism, or "
            "adverse reaction risk."
        ),
        "interpretation": (
            "Drug Response variants are relevant to pharmacogenomics and "
            "personalised medicine. They inform drug selection, dosing, or "
            "avoidance based on the patient's genotype. Clinical implementation "
            "follows guidelines from PharmGKB, CPIC, or DPWG. These variants "
            "do not indicate disease risk per se, but affect how a patient "
            "responds to specific therapeutics."
        ),
        "evidence_basis": [
            "PharmGKB clinical annotations with level of evidence",
            "CPIC or DPWG dosing guidelines referencing this variant",
            "Pharmacokinetic studies demonstrating altered drug metabolism",
            "FDA drug label pharmacogenomic information",
            "ClinVar classification as drug_response with supporting submissions",
        ],
        "acmg_criteria": None,
        "odds_ratio": None,
        "beta_value": None,
        "p_value": None,
        "population_frequency_note": None,
    },
}

# Fallback for unknown or missing classifications
_SIGNIFICANCE_FALLBACK = {
    "definition": "This classification is not recognized in the standard ClinVar terminology.",
    "interpretation": (
        "The variant's clinical significance could not be mapped to a standard "
        "category. This may indicate a novel or non-standard annotation. "
        "Manual review is recommended."
    ),
    "evidence_basis": [
        "Classification source could not be determined",
    ],
    "acmg_criteria": None,
    "odds_ratio": None,
    "beta_value": None,
    "p_value": None,
    "population_frequency_note": None,
}


def get_significance_explanation(significance: str) -> dict:
    """Return the explanation dict for a clinical significance, with fallback."""
    if not significance:
        return _SIGNIFICANCE_FALLBACK
    # Normalise common display forms to dict keys
    key = significance.strip()
    # Handle space-separated forms like "Likely pathogenic"
    normalised = key.replace(" ", "_")
    return SIGNIFICANCE_EXPLANATIONS.get(
        normalised,
        SIGNIFICANCE_EXPLANATIONS.get(key, _SIGNIFICANCE_FALLBACK),
    )


def generate_frequency_comparison(pop_freqs: dict) -> dict | None:
    """
    Calculates absolute and fold differences specifically between
    South Asian (SAS) and European (EUR) populations.
    """
    if not pop_freqs or not isinstance(pop_freqs, dict) or pop_freqs.get("available") is False:
        return None
        
    sas = pop_freqs.get("south_asian")
    eur = pop_freqs.get("european")
    
    if sas is None or eur is None:
        return None
        
    abs_diff = abs(sas - eur)
    
    # Avoid division by zero
    if min(sas, eur) == 0:
        if max(sas, eur) == 0:
            fold = 1.0
        else:
            fold = float('inf')
    else:
        fold = max(sas, eur) / min(sas, eur)
        
    return {
        "south_asian": sas,
        "european": eur,
        "absolute_difference": abs_diff,
        "fold_difference": fold,
        "higher_in": "South Asian" if sas > eur else ("European" if eur > sas else "Equal")
    }

def generate_frequency_interpretation(pop_freqs: dict) -> str:
    """
    Generates a natural language interpretation comparing South Asian and European frequencies.
    """
    comparison = generate_frequency_comparison(pop_freqs)
    
    if not comparison:
        if pop_freqs and pop_freqs.get("global") is not None:
            return f"This variant has a global allele frequency of {pop_freqs['global']:.4f}, but specific SAS/EUR ethnic breakdowns are not available."
        return "Detailed population breakdowns are not available for this variant."
        
    sas = comparison["south_asian"]
    eur = comparison["european"]
    
    if sas == eur:
        return f"This variant appears at an identical frequency ({sas:.4f}) in both South Asian and European populations."
        
    higher_name = comparison["higher_in"]
    lower_name = "European" if higher_name == "South Asian" else "South Asian"
    
    fold = comparison["fold_difference"]
    if fold == float('inf'):
        fold_text = "infinitely"
    else:
        fold_text = f"{fold:.1f}×"
        
    return (
        f"This variant exhibits a notable prevalence difference between the target populations. "
        f"It is {fold_text} more frequent in {higher_name} populations compared to {lower_name} populations, "
        f"with an absolute frequency difference of {comparison['absolute_difference']:.4f}."
    )

def generate_gwas_interpretation(assoc: dict) -> dict:
    """
    Generates an interpretation dictionary for a GWAS association based on
    P-value, Odds Ratio, and Beta.
    """
    if not assoc:
        return None
        
    p_val = assoc.get("p_value")
    or_val = assoc.get("odds_ratio")
    beta = assoc.get("beta")
    trait = assoc.get("trait", "this trait")
    
    interp = {
        "association_strength": "Unknown",
        "odds_ratio_interpretation": None,
        "beta_interpretation": None,
        "narrative": ""
    }
    
    # 1. Evaluate P-value strength
    if p_val is not None:
        if p_val < 5e-8:
            interp["association_strength"] = "Genome-wide significant"
        elif p_val < 1e-5:
            interp["association_strength"] = "Strong evidence"
        elif p_val < 1e-3:
            interp["association_strength"] = "Suggestive evidence"
        else:
            interp["association_strength"] = "Limited evidence"
            
    # 2. Evaluate Odds Ratio
    if or_val is not None:
        if or_val < 1.0:
            interp["odds_ratio_interpretation"] = "Protective"
        elif or_val < 1.5:
            interp["odds_ratio_interpretation"] = "Small increased risk"
        elif or_val < 3.0:
            interp["odds_ratio_interpretation"] = "Moderate increased risk"
        else:
            interp["odds_ratio_interpretation"] = "Strong increased risk"
            
    # 3. Evaluate Beta
    if beta is not None:
        if beta > 0:
            interp["beta_interpretation"] = "Positive association"
        elif beta < 0:
            interp["beta_interpretation"] = "Negative association"
        else:
            interp["beta_interpretation"] = "Neutral association"
            
    # 4. Construct narrative
    strength = interp["association_strength"].lower()
    if "evidence" not in strength:
        strength += " evidence"
    narrative = f"There is {strength} of an association with {trait}"
    if or_val is not None:
        narrative += f", interpreted as a {interp['odds_ratio_interpretation'].lower()} (OR: {or_val:.2f})."
    elif beta is not None:
        narrative += f", interpreted as a {interp['beta_interpretation'].lower()} (\u03b2: {beta})."
    else:
        narrative += "."
        
    interp["narrative"] = narrative
    return interp

# ═══════════════════════════════════════════════════════════════════════════════
# PubMed Literature Logic
# ═══════════════════════════════════════════════════════════════════════════════

def generate_literature_score(pubmed_data: dict) -> dict:
    """
    Computes evidence level strictly based on the total number of PubMed 
    publications returned via esearch.count.
    """
    if not pubmed_data or not pubmed_data.get("available"):
        return {"paper_count": 0, "evidence_level": "No Evidence", "score": 0}
        
    count = pubmed_data.get("paper_count", 0)
    
    if count == 0:
        level = "No Evidence"
        score = 0
    elif 1 <= count <= 10:
        level = "Limited Evidence"
        score = 1
    elif 11 <= count <= 100:
        level = "Emerging Evidence"
        score = 2
    elif 101 <= count <= 1000:
        level = "Moderate Evidence"
        score = 3
    else:
        level = "Strong Evidence"
        score = 4
        
    return {
        "paper_count": count,
        "evidence_level": level,
        "score": score
    }

def generate_literature_interpretation(pubmed_data: dict) -> str:
    """
    Provides a standardized text interpretation of the literature volume
    without making any biological conclusions.
    """
    score_data = generate_literature_score(pubmed_data)
    level = score_data["evidence_level"]
    
    if level == "No Evidence":
        return "No relevant publications were identified."
    elif level == "Limited Evidence":
        return "Limited literature is currently available."
    elif level == "Emerging Evidence":
        return "Emerging literature exists, suggesting ongoing research into this variant."
    elif level == "Moderate Evidence":
        return "A moderate volume of literature is available, indicating established research interest."
    else:
        return "Extensive literature exists supporting investigation of this variant."

# ═══════════════════════════════════════════════════════════════════════════════
# INTERPRETATION SUMMARY ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

_DISCLAIMER = (
    "This interpretation is intended for research and educational purposes "
    "and should not be considered clinical advice."
)

# Impact-level → terse functional description (used inside the narrative)
_IMPACT_PHRASES = {
    "HIGH": (
        "predicted to substantially alter or abolish biological function — "
        "consistent with loss-of-function"
    ),
    "MODERATE": (
        "predicted to alter protein structure or activity without completely "
        "disrupting the gene product"
    ),
    "LOW": (
        "predicted to have minimal or no effect on protein function"
    ),
    "MODIFIER": (
        "located outside protein-coding sequence with uncertain functional relevance"
    ),
}

# Evidence-strength → qualifier phrase
_EVIDENCE_PHRASES = {
    "High":     "supported by high-confidence evidence",
    "Moderate": "supported by moderate-confidence evidence",
    "Low":      "supported by limited evidence",
}


def generate_gene_relevance_score(pubmed_count: int, disease_count: int, has_clinvar: bool, has_gwas: bool) -> str:
    """
    Generate a Research Relevance score based on available biological and clinical evidence.
    """
    score = 0
    if pubmed_count > 1000:
        score += 3
    elif pubmed_count > 100:
        score += 2
    elif pubmed_count > 10:
        score += 1
        
    if disease_count > 5:
        score += 2
    elif disease_count > 0:
        score += 1
        
    if has_clinvar:
        score += 2
        
    if has_gwas:
        score += 1
        
    if score >= 6:
        return "Very High"
    elif score >= 4:
        return "High"
    elif score >= 2:
        return "Moderate"
    else:
        return "Low"


def generate_interpretation_summary(
    gene: str,
    impact: str,
    clinical_significance: str,
    condition: str | None,
    evidence_score: int,
    evidence_strength: str,
    variant_type: str = "SNP",
    clinvar_evidence: dict | None = None,
) -> dict:
    """
    Synthesise a human-readable scientific interpretation narrative that
    explains the relationship between variant impact, clinical significance,
    gene context, and available evidence.

    Parameters
    ----------
    gene                 : Gene symbol (e.g. "BRCA1").  May be None for intergenic.
    impact               : Impact level string — "HIGH" | "MODERATE" | "LOW" | "MODIFIER"
    clinical_significance: ClinVar classification string (or None / "Not assessed").
    condition            : Associated disease / condition string (may be None).
    evidence_score       : Integer 0–100 from the evidence framework.
    evidence_strength    : "High" | "Moderate" | "Low" from the evidence framework.
    variant_type         : "SNP" | "CNV" | "Indel" (for phrasing variation).

    Returns
    -------
    dict with a single key:
        interpretation_summary : str — multi-sentence narrative + disclaimer.

    Design principles
    -----------------
    - Generates classification-specific paragraphs: each ClinVar family has a
      distinct narrative template.
    - Integrates impact × significance cross-commentary when the two convey
      different signals (e.g. HIGH impact but Risk_factor significance).
    - Appends the mandatory research disclaimer as a separate final sentence.
    - Does NOT modify any classification, impact score, or evidence score.
    - Reuses phrasing from IMPACT_EXPLANATIONS and SIGNIFICANCE_EXPLANATIONS
      dictionaries so the narrative is always semantically consistent with the
      expandable explanation cards.
    """
    # ── Normalise inputs ───────────────────────────────────────────────────────
    gene_name   = gene.upper() if gene else "an unknown gene"
    impact_norm = (impact or "MODIFIER").upper()
    sig_raw     = (clinical_significance or "").strip()
    sig_key     = sig_raw.replace(" ", "_")
    condition   = condition or None
    ev_strength = evidence_strength or "Low"
    ev_score    = evidence_score or 0
    vtype       = variant_type or "variant"

    impact_phrase  = _IMPACT_PHRASES.get(impact_norm, _IMPACT_PHRASES["MODIFIER"])
    evidence_phrase = _EVIDENCE_PHRASES.get(ev_strength, "with available evidence")

    # Condition clause
    condition_clause = f" associated with {condition}" if condition else ""

    # ── Opening sentence — always present ─────────────────────────────────────
    opening = (
        f"This {vtype.lower()} in {gene_name} is {impact_phrase}. "
        f"The overall evidence base yields a score of {ev_score}/100, "
        f"{evidence_phrase}."
    )

    # ── Classification-specific body paragraph ─────────────────────────────────
    if sig_key in ("Pathogenic", "Likely_pathogenic"):
        certainty = "current evidence strongly supports" if sig_key == "Pathogenic" \
                    else "current evidence suggests — though does not definitively confirm —"
        body = (
            f"{certainty.capitalize()} a disease-causing role for this variant"
            f"{condition_clause}. "
        )
        if sig_key == "Likely_pathogenic":
            body += (
                "A Likely Pathogenic classification indicates that the evidence "
                "is compelling but not yet definitive; reclassification may occur "
                "as additional functional or segregation data become available. "
            )
        if impact_norm == "HIGH":
            body += (
                "The HIGH functional impact prediction is consistent with the "
                "pathogenic classification: high-impact variants commonly disrupt "
                "gene function through nonsense-mediated decay, frameshifts, or "
                "critical splice-site disruption."
            )
        elif impact_norm == "MODERATE":
            body += (
                "The MODERATE functional impact prediction is consistent with a "
                "protein-altering missense or in-frame change that may impair, "
                "but does not necessarily abolish, normal protein function."
            )
        else:
            body += (
                f"Note that the {impact_norm} functional impact prediction may appear "
                "discordant with the pathogenic classification; this can occur when "
                "the variant exerts its effect through a regulatory or dosage-sensitive "
                "mechanism not captured by standard consequence predictions."
            )

    elif sig_key == "Risk_factor":
        body = (
            f"This variant is classified as a risk factor{condition_clause}, "
            "meaning it statistically increases disease susceptibility without "
            "being necessary or sufficient to cause disease independently. "
        )
        if impact_norm == "HIGH":
            body += (
                "Although the variant is predicted to have HIGH functional impact, "
                "its clinical classification as a Risk Factor indicates that carrying "
                "the variant does not deterministically lead to disease — penetrance "
                "and expressivity are modulated by additional genetic and environmental factors."
            )
        elif impact_norm in ("LOW", "MODIFIER"):
            body += (
                f"The {impact_norm} predicted functional impact is consistent with a "
                "risk-modifying role: the variant likely exerts a subtle quantitative "
                "effect on gene expression or protein dosage rather than an overt "
                "loss-of-function."
            )
        else:
            body += (
                "Risk stratification for this variant should be interpreted in the "
                "context of polygenic risk scores, family history, and relevant "
                "environmental exposures."
            )

    elif sig_key in ("Benign", "Likely_benign"):
        certainty = "established" if sig_key == "Benign" else "likely"
        body = (
            f"Available evidence does not support a disease-causing role for this "
            f"variant; it is {certainty} benign"
            f"{condition_clause}. "
        )
        if sig_key == "Likely_benign":
            body += (
                "A Likely Benign classification carries a small residual uncertainty "
                "and warrants periodic re-evaluation if the clinical phenotype is "
                "highly atypical. "
            )
        body += (
            f"The {impact_norm} predicted functional impact is consistent with "
            f"a common polymorphism or neutral coding change that does not "
            f"appreciably alter gene function."
        )

    elif sig_key == "Uncertain_significance":
        body = (
            "Evidence is currently insufficient to classify this variant as "
            "either disease-causing or benign (Variant of Uncertain Significance, "
            "VUS). "
            "This classification should NOT be used for clinical decision-making "
            "without additional supporting data. "
        )
        if impact_norm == "HIGH":
            body += (
                "The HIGH predicted functional impact increases the prior probability "
                "of pathogenicity and may support prioritisation for functional "
                "follow-up studies, but alone is insufficient to reclassify the variant."
            )
        else:
            body += (
                f"The {impact_norm} predicted functional impact does not provide "
                "strong evidence in either direction. Periodic reinterpretation through "
                "ClinVar or the originating laboratory is strongly recommended."
            )

    elif sig_key in ("Drug_response", "drug_response"):
        body = (
            f"This variant in {gene_name} is associated with altered pharmacological "
            f"response{condition_clause}, rather than with disease causation per se. "
            "Drug response variants are relevant to personalised medicine: they may "
            "affect drug efficacy, metabolism rate, or adverse reaction risk and are "
            "managed according to pharmacogenomics guidelines (CPIC / PharmGKB). "
        )
        if impact_norm in ("MODERATE", "HIGH"):
            body += (
                f"The {impact_norm} predicted functional impact is consistent with "
                "a meaningful change in the encoded protein's enzymatic or receptor "
                "activity that underlies the altered drug response phenotype."
            )
        else:
            body += (
                "The variant may exert its pharmacogenomic effect through subtle "
                "changes in gene expression or protein stability rather than "
                "a dramatic loss-of-function."
            )

    else:
        # Fallback for "Not assessed", "Not in ClinVar", or any unrecognised classification
        body = (
            f"No clinical classification is currently available for this variant "
            f"in {gene_name}. "
            f"The {impact_norm} predicted functional impact provides a preliminary "
            f"indication of potential biological consequence, but clinical significance "
            f"cannot be determined without formal variant interpretation. "
            f"Submission to ClinVar and consultation with a certified clinical "
            f"genomics laboratory is recommended."
        )

    # ── ClinVar Evidence Engine clause ─────────────────────────────────────────
    clinvar_clause = ""
    if clinvar_evidence:
        rev_status = clinvar_evidence.get("review_status", "no assertion")
        conf_level = clinvar_evidence.get("confidence_level", "Low")
        cons_score = clinvar_evidence.get("consensus_score")
        
        cv_text = []
        if "expert panel" in rev_status.lower() or "practice guideline" in rev_status.lower():
            cv_text.append("is supported by expert panel review")
        else:
            cv_text.append(f"has a {conf_level.lower()} confidence level ({rev_status})")
            
        if cons_score is not None:
            if cons_score >= 80:
                cv_text.append("demonstrates high consensus among current clinical submissions")
            elif cons_score >= 50:
                cv_text.append("demonstrates moderate consensus among submissions")
            else:
                cv_text.append("has highly conflicting clinical submissions")
                
        if cv_text:
            clinvar_clause = f" The ClinVar classification {' and '.join(cv_text)}."

    # ── Final assembly ─────────────────────────────────────────────────────────
    summary = f"{opening}{body}{clinvar_clause} {_DISCLAIMER}"

    return {"interpretation_summary": summary}


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
        "found": True,
        "error": None,
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
        "evidence": build_evidence_object(None, gene_record["gene"] if gene_record else None),
        "population_frequencies": None,
        "frequency_source": "Not applicable for unannotated coordinate lookup",
        "frequency_interpretation": "Detailed population breakdowns are not available for unannotated coordinates.",
        "gwas_associations": None,
        "gwas_top_hit": None,
        "gwas_interpretation": None,
        "impact_explanation": get_impact_explanation(impact.get("impact_level")),
        "significance_explanation": get_significance_explanation(None),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 3b. SNP ANNOTATION FROM rsID
# ═══════════════════════════════════════════════════════════════════════════════

def annotate_snp_from_rsid(rsid: str) -> dict:
    """
    Full annotation pipeline entered via rsID rather than coordinates.

    Resolves chromosome / position / ref / alt from the mock dbSNP dataset,
    then runs the identical interpretation chain as annotate_snp() and enriches
    the result with ClinVar clinical significance + condition.

    Returns the same dict shape as annotate_snp() plus:
        rsid               : the queried rsID
        source             : "rsid_lookup"
        hgvs               : HGVS notation from dbSNP
        protein_change     : protein-level change (if missense)
        allele_frequency   : global allele frequency
        clinical_significance : ClinVar classification (may be None)
        condition          : associated disease/condition (may be None)
        review_status      : ClinVar review status (may be None)
    """
    snp_record     = get_dbsnp_record(rsid)
    clinvar_record = snp_record.get("clinvar") if snp_record else None

    if not snp_record:
        return {
            "found": False,
            "error": f"rsID '{rsid}' was not found in the local dataset. "
                     "In production this would query NCBI E-utilities.",
            "rsid": rsid,
        }

    chromosome  = str(snp_record.get("chromosome", ""))
    position    = snp_record.get("position", 0)
    ref         = snp_record.get("ref", "A")
    alt         = snp_record.get("alt", "G")
    consequence = snp_record.get("consequence")

    # Reuse existing pipeline functions — nothing is duplicated
    classification = classify_variant(chromosome, position, ref, alt)
    raw_gene = snp_record.get("gene", "")
    gene_record = get_gene_info(raw_gene)
    if not gene_record and raw_gene:
        # Dynamic fallback if the gene is not in gene_coordinates.json
        gene_record = {
            "gene": raw_gene,
            "omim": "Unknown",
            "pathway": "Unknown",
            "chromosome": chromosome,
            "start": position,
            "end": position,
            "strand": "+"
        }
        
    impact      = predict_functional_impact(ref, alt, consequence)
    gene_info   = gene_record          # same object; alias for clarity
    chrom_info  = get_chromosome_info(chromosome)

    clinical_significance = (
        clinvar_record.get("clinical_significance") if clinvar_record else None
    )

    gwas_assocs = snp_record.get("gwas_associations")

    pubmed_data = _fetch_live_pubmed(rsid, gene_record.get("gene") if gene_record else "")
    pubmed_score = generate_literature_score(pubmed_data)
    pubmed_interp = generate_literature_interpretation(pubmed_data)
    pubmed_data.update(pubmed_score)
    pubmed_data["interpretation"] = pubmed_interp

    # Version 2.2: Gene Context & Biological Annotation Engine
    gene_symbol_clean = gene_record.get("gene") if gene_record else None
    if gene_symbol_clean:
        gene_context = _fetch_live_gene_context(gene_symbol_clean)
        diseases = _fetch_gene_diseases(gene_symbol_clean)
        pathways = _fetch_gene_pathways(gene_symbol_clean)
    else:
        gene_context = {"available": False}
        diseases = {"available": False, "diseases": []}
        pathways = {"available": False, "pathways": []}
        
    research_relevance = generate_gene_relevance_score(
        pubmed_count=pubmed_data.get("paper_count", 0),
        disease_count=len(diseases.get("diseases", [])),
        has_clinvar=bool(clinical_significance is not None and clinical_significance != "Unknown"),
        has_gwas=bool(gwas_assocs)
    )

    return {
        "found": True,
        "error": None,
        "input": {
            "chromosome": chromosome,
            "position":   position,
            "ref":        ref.upper(),
            "alt":        alt.upper(),
        },
        "rsid":          rsid,
        "classification": classification,
        "gene":          gene_record,
        "gene_info":     gene_info,
        "impact":        impact,
        "chromosome_info": chrom_info,
        "source":        "rsid_lookup",
        # rsID-enriched fields
        "hgvs":                snp_record.get("hgvs"),
        "protein_change":      snp_record.get("protein_change"),
        "allele_frequency":    snp_record.get("allele_frequency"),
        "population_frequencies": snp_record.get("population_frequencies"),
        "frequency_source":    snp_record.get("frequency_source"),
        "frequency_interpretation": generate_frequency_interpretation(snp_record.get("population_frequencies")),
        "frequency_comparison": generate_frequency_comparison(snp_record.get("population_frequencies")),
        "gwas_associations":   gwas_assocs,
        "gwas_top_hit":        gwas_assocs[0] if gwas_assocs else None,
        "gwas_interpretation": generate_gwas_interpretation(gwas_assocs[0] if gwas_assocs else None),
        "clinical_significance": clinical_significance,
        "condition":           clinvar_record.get("condition") if clinvar_record else None,
        "review_status":       clinvar_record.get("review_status") if clinvar_record else None,
        "evidence": build_evidence_object(rsid, snp_record.get("gene")),
        "impact_explanation":       get_impact_explanation(impact.get("impact_level")),
        "significance_explanation": get_significance_explanation(clinical_significance),
        "interpretation_summary":   generate_interpretation_summary(
            gene=gene_record.get("gene") if gene_record else None,
            impact=impact.get("impact_level"),
            clinical_significance=clinical_significance,
            condition=clinvar_record.get("condition") if clinvar_record else None,
            evidence_score=build_evidence_object(rsid, snp_record.get("gene")).get("evidence_score", 0),
            evidence_strength=build_evidence_object(rsid, snp_record.get("gene")).get("evidence_strength", "Low"),
            variant_type="SNP",
            clinvar_evidence=clinvar_record
        ).get("interpretation_summary"),
        "clinvar": clinvar_record,
        "pubmed": pubmed_data,
        "gene_context": gene_context,
        "diseases": diseases,
        "pathways": pathways,
        "research_relevance": research_relevance
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
        "found": True,
        "error": None,
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
        "evidence": build_evidence_object(None, gene_record["gene"] if gene_record else None),
        "population_frequencies": None,
        "frequency_source": "Not applicable for structural variants",
        "frequency_interpretation": "Detailed population breakdowns are not available for structural variants.",
        "gwas_associations": None,
        "gwas_top_hit": None,
        "gwas_interpretation": None,
        "impact_explanation": get_impact_explanation(dosage_class),
        "significance_explanation": get_significance_explanation(None),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 4b. CNV ANALYSIS FROM rsID
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_cnv_from_rsid(rsid: str, cnv_type: str,
                           copy_number: int = None) -> dict:
    """
    Resolve an rsID to genomic coordinates, then run CNV analysis.

    Strategy:
      1. Look up the rsID in mock dbSNP → get the gene symbol.
      2. Look up the gene in gene_coordinates.json → get its full genomic span.
      3. Use that span as the CNV region (real gene boundaries, not a ±50 kb hack).
      4. Fall back to position ± 50,000 bp only if the gene is unknown.
      5. Delegate entirely to analyze_cnv() — nothing duplicated.

    The returned dict is identical to analyze_cnv() output with two additions:
        input.rsid   : the queried rsID
        input.source : "rsid_lookup"
    """
    snp_record = get_dbsnp_record(rsid)
    if not snp_record:
        return {
            "found": False,
            "error": f"rsID '{rsid}' was not found in the local dataset.",
            "rsid": rsid,
        }

    chromosome   = str(snp_record.get("chromosome", ""))
    position     = snp_record.get("position", 0)
    gene_symbol  = snp_record.get("gene", "")

    # Use real gene boundaries when available
    gene_info = get_gene_info(gene_symbol) if gene_symbol else None
    if gene_info:
        start = gene_info["start"]
        end   = gene_info["end"]
    else:
        start = max(1, position - 50_000)
        end   = position + 50_000

    result = analyze_cnv(chromosome, start, end, cnv_type, copy_number)
    result["input"]["rsid"]   = rsid
    result["input"]["source"] = "rsid_lookup"
    result["clinvar"]         = snp_record.get("clinvar")
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# 4c. CNV ANALYSIS FROM GENE SYMBOL
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_cnv_from_gene(gene_symbol: str, cnv_type: str,
                           copy_number: int = None) -> dict:
    """
    Resolve a gene symbol to its genomic coordinates, then run CNV analysis.

    Uses the full gene span (start → end from gene_coordinates.json) so the
    CNV region precisely covers the targeted gene body.

    The returned dict is identical to analyze_cnv() output with two additions:
        input.gene_symbol : the queried gene symbol
        input.source      : "gene_lookup"
    """
    gene_info = get_gene_info(gene_symbol)
    if not gene_info:
        return {
            "found": False,
            "error": f"Gene '{gene_symbol}' was not found in the local dataset.",
        }

    chromosome = gene_info["chromosome"]
    start      = gene_info["start"]
    end        = gene_info["end"]

    result = analyze_cnv(chromosome, start, end, cnv_type, copy_number)
    result["input"]["gene_symbol"] = gene_symbol
    result["input"]["source"]      = "gene_lookup"
    return result


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
        clinvar_record = snp_record.get("clinvar") if snp_record else None

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
                "error": None,
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
                "evidence": build_evidence_object(rsid, snp_record.get("gene")),
                "impact_explanation": get_impact_explanation(impact.get("impact_level")),
                "significance_explanation": get_significance_explanation(
                    clinvar_record.get("clinical_significance") if clinvar_record else None
                ),
            })
        else:
            results.append({
                "rsid": rsid,
                "found": False,
                "error": f"rsID '{rsid}' not found.",
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
                "evidence": build_evidence_object(rsid, None),
                "impact_explanation": get_impact_explanation("MODIFIER"),
                "significance_explanation": get_significance_explanation("Not in database"),
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
    base_report = annotate_snp_from_rsid(rsid)
    if not base_report.get("found"):
        return {"found": False, "rsid": rsid}

    # Re-fetch snp_record to provide the 'snp' key expected by report.html
    base_report["snp"] = get_dbsnp_record(rsid)
    
    snp_record     = base_report.get("snp", {})
    clinvar_record = base_report.get("clinvar")
    impact         = base_report.get("impact", {})
    classification = base_report.get("classification", {})
    gene_info      = base_report.get("gene_info", {})

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

    interpretation = {
        "impact": impact.get("impact_level", "UNKNOWN"),
        "impact_meaning": impact.get("description", "Meaning unknown."),
        "clinical_significance": significance,
        "significance_meaning": (
            "Existing evidence suggests disease association." if significance in ["Pathogenic", "Likely_pathogenic"]
            else "Not enough evidence for pathogenicity." if significance == "Uncertain_significance"
            else "Considered non-pathogenic." if significance in ["Benign", "Likely_benign"]
            else "Variant is a recognized risk factor for disease but not deterministically pathogenic." if significance == "Risk_factor"
            else "Variant affects pharmacological response to specific drugs." if significance == "drug_response"
            else "Not assessed."
        ),
        "general_summary": (
            f"Variant {rsid} ({hgvs}) introduces a {classification.get('subtype', '').lower()} "
            f"in {gene_name}, resulting in a {consequence.replace('_', ' ')} "
            f"({protein_change}). "
            f"The global alternate allele frequency is {freq_str}. "
            f"{'It participates in the ' + gene_info.get('pathway', 'unknown') + ' pathway.' if gene_info else ''}"
        )
    }

    base_report["interpretation"] = interpretation
    base_report["sig_colour"] = sig_colour
    
    return base_report
