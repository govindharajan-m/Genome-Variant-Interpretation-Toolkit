# Genome Variation Analysis Platform (GenomeVAP)

**Version:** 3.0
**Author:** B.Tech Industrial Biotechnology Portfolio Project

GenomeVAP is an advanced, terminal-aesthetic web platform built for genomic variant interpretation. It mimics a 1990s Unix lab workstation while providing modern, complex analytical pipelines for bioinformatics analysis. 

## Overview
GenomeVAP transforms raw genomic data (rsIDs, chromosomal coordinates) into comprehensive, biologically contextualized reports. The platform leverages a deterministic scoring engine to assess variant priority, discovery potential, and research relevance without relying on stochastic models.

## Features
- **Single Variant Annotation (SNP/Indel):** Deep-dive analysis of specific loci.
- **Batch rsID Analysis:** High-throughput triage of multiple variants.
- **Copy Number Variant (CNV) Analysis:** Structural variant interpretation.
- **Cohort Analysis:** Population-level variant prioritization.
- **Comparative Analysis:** Head-to-head evaluation of candidate variants.
- **Disease Panel Designer:** Custom gene panel recommendations based on discovery potential.
- **Pathway Analysis:** Systems biology mapping of variants to functional pathways.
- **Client-Side Export:** Secure, zero-latency CSV reporting.

## Architecture
GenomeVAP is a Python/Flask monolith.
- **Backend:** Flask routes handle request parsing and delegation.
- **Engine:** `variant_engine.py` executes all mathematical and biological scoring heuristics.
- **Database:** `db_handler.py` manages hydration from localized mock datasets.
- **Frontend:** Jinja2 templates, CSS custom properties, and vanilla JS (augmented with Bootstrap 5 for responsiveness).

## Installation
```bash
# Clone the repository
git clone https://github.com/yourusername/GenomeVAP.git
cd GenomeVAP

# Install dependencies
pip install -r requirements.txt

# Run the application
python app.py
```
Access the platform at `http://localhost:5000`.

## Usage
1. Navigate to the desired tool via the dashboard.
2. Enter a variant (e.g., `rs334`, `rs429358`, `rs1042522`).
3. View the generated report or export the summary via CSV.

## Validation
All scoring engines undergo strict deterministic validation testing. The V3 release ensures 100% computational reproducibility across the Clinical, Research, and Discovery axes. See `validate_v29.py` for testing frameworks.

## Limitations
**RESEARCH USE ONLY.** GenomeVAP utilizes static mock data for demonstration purposes and simplified ACMG heuristics. See documentation for details.

## Future Work
- Integration with live external APIs (Ensembl REST, ClinVar).
- Implementation of rigorous ACMG/AMP 2015 boolean logic.
- VCF file parsing for cohort data upload.
