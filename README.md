# 🧬 GenomeVAP — Genome Variation Analysis Platform

A **portfolio-grade bioinformatics web application** built with Python (Flask) + HTML5/CSS/JS.

Simulates a simplified version of real-world tools:
- Ensembl Variant Effect Predictor (VEP)
- ANNOVAR
- ClinVar annotation pipelines

---

## 📸 Features

| Module | Description |
|--------|-------------|
| **SNP Analysis** | Enter Chr:Pos:Ref:Alt → gene mapping, consequence prediction (missense/nonsense/synonymous), SIFT & PolyPhen-2 |
| **CNV Analysis** | Deletion/Duplication region → dosage effect, haploinsufficiency, gene affected |
| **Batch rsID** | Up to 200 rsIDs → full table with gene, consequence, ClinVar significance, CSV export |
| **Variant Reports** | Full per-rsID report: dbSNP + ClinVar + gene info + scientific interpretation |

---

## 🗂 Project Structure

```
genome_variant_platform/
├── app.py                  # Flask routes & server
├── variant_engine.py       # Bioinformatics analysis logic
├── db_handler.py           # Mock dbSNP / ClinVar / Ensembl data handler
├── requirements.txt
├── data/
│   ├── dbsnp_mock.json     # 15 real rsIDs with coordinates & consequences
│   ├── clinvar_mock.json   # Clinical significance for all 15 rsIDs
│   └── gene_coordinates.json  # 14 genes with genomic coordinates
├── templates/
│   ├── base.html           # Layout with navbar & footer
│   ├── index.html          # Home / landing page
│   ├── single_variant.html # SNP analysis page
│   ├── cnv_analysis.html   # CNV analysis page
│   ├── batch_analysis.html # Batch rsID page
│   ├── report.html         # Full variant report page
│   └── 404.html
└── static/
    ├── css/style.css       # Dark bioinformatics UI
    └── js/main.js          # Shared JS utilities
```

---

## 🚀 Quick Start

### 1. Clone / Download the project

```bash
cd genome_variant_platform
```

### 2. Create a virtual environment (recommended)

```bash
python -m venv venv
source venv/bin/activate          # Linux / macOS
venv\Scripts\activate             # Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Run the server

```bash
python app.py
```

Open **http://localhost:5000** in your browser.

---

## 🧪 Demo Variants (included in mock dataset)

| rsID | Gene | Consequence | Clinical Significance |
|------|------|-------------|----------------------|
| rs334 | HBB | Missense (p.Glu7Val) | **Pathogenic** — Sickle Cell Disease |
| rs28897696 | CFTR | Stop gained | **Pathogenic** — Cystic Fibrosis |
| rs1800562 | HFE | Missense | **Pathogenic** — Hereditary Hemochromatosis |
| rs113488022 | BRAF | Missense (V600E) | **Pathogenic** — Melanoma |
| rs28934578 | BRCA1 | Frameshift | **Pathogenic** — HBOC syndrome |
| rs80357906 | BRCA2 | Frameshift | **Pathogenic** — HBOC syndrome |
| rs429358 | APOE | Missense | Risk factor — Alzheimer's |
| rs7412 | APOE | Missense | Likely pathogenic |
| rs1042522 | TP53 | Missense (Pro72Arg) | Benign |
| rs762551 | CYP1A2 | Synonymous | Benign — Caffeine metabolism |
| rs9939609 | FTO | Intron | Risk factor — Obesity |
| rs1805007 | MC1R | Missense | Likely pathogenic — Melanoma risk |
| rs1799971 | OPRM1 | Missense | Benign — Opioid sensitivity |
| rs699 | AGT | Missense | Benign — Hypertension susceptibility |
| rs2230199 | C3 | Missense | VUS — AMD |

---

## 🔌 API Endpoints

### Single SNP Annotation
```http
POST /api/analyze-snp
Content-Type: application/json

{
  "chromosome": "17",
  "position": 7676154,
  "ref": "C",
  "alt": "G"
}
```

### CNV Analysis
```http
POST /api/analyze-cnv
Content-Type: application/json

{
  "chromosome": "17",
  "start": 43044295,
  "end": 43125483,
  "cnv_type": "Deletion",
  "copy_number": 1
}
```

### Batch rsID Analysis
```http
POST /api/batch
Content-Type: application/json

{
  "rsids": ["rs334", "rs7412", "rs429358"]
}
```

### Full Variant Report
```http
GET /api/report/rs334
```

### CSV Download
```http
POST /download/csv
Content-Type: application/json

{"results": [...]}   # same array returned by /api/batch
```

---

## 🧬 Bioinformatics Concepts Demonstrated

- **SNP Classification**: Transition vs. Transversion based on nucleotide change
- **Functional Consequence Prediction**: Using codon table simulation (synonymous / missense / nonsense)
- **CNV Dosage Analysis**: Haploinsufficiency and triplosensitivity prediction
- **Clinical Significance**: Pathogenic / Likely Pathogenic / VUS / Benign — per ACMG guidelines
- **HGVS Nomenclature**: Human Genome Variation Society standard notation
- **Gene Mapping**: Coordinate-based overlap detection (simplified interval scan)
- **dbSNP / ClinVar / Ensembl**: Data model simulation with extensible API hooks

---

## 🔧 Extending to Real APIs

Each function in `db_handler.py` has a comment with the real API equivalent:

```python
# Real-world equivalent:
#   GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
#       ?db=snp&id=<rsid>&rettype=json
```

To connect to real databases, replace the local JSON lookups with the corresponding HTTP calls using the `requests` library.

---

## 📚 Technologies

- **Backend**: Python 3.11+, Flask 3.x
- **Frontend**: HTML5, CSS3 (Flexbox/Grid), Vanilla JavaScript
- **Fonts**: Inter + JetBrains Mono (Google Fonts)
- **Data**: Local JSON (inspired by dbSNP, ClinVar, Ensembl GRCh38)
- **Reference Genome**: GRCh38 / hg38

---

*Built as a B.Tech Industrial Biotechnology portfolio project demonstrating bioinformatics pipeline concepts.*
