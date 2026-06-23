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

## Screenshots

### Single Variant Dashboard
![Single Variant Dashboard](C:/Users/mgovi/.gemini/antigravity-ide/brain/5ba47b3c-9e65-4ee7-9310-47ceacac79b5/single_variant_dashboard_1782120685917.png)

### Clinical Evidence Report
![Clinical Evidence Report View](C:/Users/mgovi/.gemini/antigravity-ide/brain/5ba47b3c-9e65-4ee7-9310-47ceacac79b5/report_view_1782120701687.png)

### Disease Panel Designer
![Disease Panel Designer](C:/Users/mgovi/.gemini/antigravity-ide/brain/5ba47b3c-9e65-4ee7-9310-47ceacac79b5/disease_panel_designer_1782120716458.png)

### Cohort Analysis
![Cohort Analysis Dashboard](C:/Users/mgovi/.gemini/antigravity-ide/brain/5ba47b3c-9e65-4ee7-9310-47ceacac79b5/cohort_analysis_1782120733700.png)

### Comparative Analysis
![Comparative Analysis Grid](C:/Users/mgovi/.gemini/antigravity-ide/brain/5ba47b3c-9e65-4ee7-9310-47ceacac79b5/comparative_analysis_1782120750014.png)

### Pathway Analysis
![Systems Biology Pathway Mapping](C:/Users/mgovi/.gemini/antigravity-ide/brain/5ba47b3c-9e65-4ee7-9310-47ceacac79b5/pathway_analysis_1782120763437.png)


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
