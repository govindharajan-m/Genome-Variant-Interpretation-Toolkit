"""
app.py — Flask Web Application Entry Point
============================================
Genome Variation Analysis Platform
Built for bioinformatics portfolio demonstration.

Architecture:
  - /                        → Home dashboard
  - /single-variant          → SNP / Indel input form
  - /cnv-analysis            → CNV input form
  - /batch-analysis          → Batch rsID input form
  - /api/analyze-snp         → POST: annotate a single SNP
  - /api/analyze-cnv         → POST: annotate a CNV
  - /api/batch               → POST: batch rsID analysis
  - /api/report/<rsid>       → GET: full report for one rsID
  - /download/csv            → POST: export batch results as CSV
"""

import csv
import io
import json

from flask import (
    Flask, render_template, request,
    jsonify, make_response, redirect, url_for, session
)
from variant_engine import (
    annotate_snp,
    analyze_cnv,
    batch_analyze_rsids,
    generate_interpretation_report,
    predict_functional_impact,
    classify_variant,
    SIGNIFICANCE_COLOURS,
)
from db_handler import get_all_known_rsids

app = Flask(__name__)
app.secret_key = "gvap_demo_secret_2024"   # change for production


# ═══════════════════════════════════════════════════════════════════════════════
# TEMPLATE CONTEXT HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

@app.context_processor
def inject_globals():
    """Inject values available to all Jinja2 templates."""
    return {
        "sig_colours": SIGNIFICANCE_COLOURS,
        "platform_name": "GenomeVAP",
        "platform_version": "1.0",
    }


# ═══════════════════════════════════════════════════════════════════════════════
# PAGE ROUTES
# ═══════════════════════════════════════════════════════════════════════════════

@app.route("/")
def home():
    """Landing page with platform overview and quick-access cards."""
    known_rsids = get_all_known_rsids()[:6]   # show first 6 as demo examples
    return render_template("index.html", demo_rsids=known_rsids)


@app.route("/single-variant")
def single_variant_page():
    """Page for manual SNP / Indel annotation by chromosomal coordinates."""
    return render_template("single_variant.html")


@app.route("/cnv-analysis")
def cnv_analysis_page():
    """Page for CNV (Copy Number Variant) analysis."""
    return render_template("cnv_analysis.html")


@app.route("/batch-analysis")
def batch_analysis_page():
    """Page for batch rsID annotation."""
    known_rsids = get_all_known_rsids()
    return render_template("batch_analysis.html", known_rsids=known_rsids)


@app.route("/report/<rsid>")
def report_page(rsid):
    """Full detailed report page for a single rsID."""
    report = generate_interpretation_report(rsid.lower())
    return render_template("report.html", report=report, rsid=rsid)


# ═══════════════════════════════════════════════════════════════════════════════
# API ENDPOINTS (JSON)
# ═══════════════════════════════════════════════════════════════════════════════

@app.route("/api/analyze-snp", methods=["POST"])
def api_analyze_snp():
    """
    Annotate a single SNP given chromosomal coordinates.

    Expected JSON body:
        { "chromosome": "17", "position": 7676154, "ref": "C", "alt": "G" }

    Returns:
        JSON with full annotation result.
    """
    data = request.get_json(force=True)

    # ── Input validation ───────────────────────────────────────────────────────
    required = ["chromosome", "position", "ref", "alt"]
    missing  = [f for f in required if not data.get(f)]
    if missing:
        return jsonify({"error": f"Missing fields: {', '.join(missing)}"}), 400

    try:
        position = int(data["position"])
    except (ValueError, TypeError):
        return jsonify({"error": "Position must be an integer"}), 400

    valid_bases = set("ACGTN")
    for field in ("ref", "alt"):
        if not set(data[field].upper()).issubset(valid_bases):
            return jsonify({"error": f"Invalid nucleotide in {field}"}), 400

    result = annotate_snp(
        chromosome=data["chromosome"].upper().replace("CHR", ""),
        position=position,
        ref=data["ref"],
        alt=data["alt"],
    )

    return jsonify(result)


@app.route("/api/analyze-cnv", methods=["POST"])
def api_analyze_cnv():
    """
    Annotate a CNV (Copy Number Variant).

    Expected JSON body:
        {
            "chromosome": "17",
            "start": 43044295,
            "end": 43125483,
            "cnv_type": "Deletion",
            "copy_number": 1
        }
    """
    data = request.get_json(force=True)

    required = ["chromosome", "start", "end", "cnv_type"]
    missing  = [f for f in required if not data.get(f)]
    if missing:
        return jsonify({"error": f"Missing fields: {', '.join(missing)}"}), 400

    try:
        start = int(data["start"])
        end   = int(data["end"])
    except (ValueError, TypeError):
        return jsonify({"error": "start and end must be integers"}), 400

    if start >= end:
        return jsonify({"error": "start must be less than end"}), 400

    cnv_type = data["cnv_type"]
    if cnv_type.lower() not in ("duplication", "deletion"):
        return jsonify({"error": "cnv_type must be 'Duplication' or 'Deletion'"}), 400

    copy_number = data.get("copy_number")
    if copy_number is not None:
        try:
            copy_number = int(copy_number)
        except ValueError:
            return jsonify({"error": "copy_number must be an integer"}), 400

    result = analyze_cnv(
        chromosome=data["chromosome"].upper().replace("CHR", ""),
        start=start,
        end=end,
        cnv_type=cnv_type,
        copy_number=copy_number,
    )

    return jsonify(result)


@app.route("/api/batch", methods=["POST"])
def api_batch():
    """
    Batch rsID annotation.

    Accepts either:
      - JSON body: { "rsids": ["rs334", "rs7412", ...] }
      - Multipart form with field "rsid_list" (text) or file upload "rsid_file"

    Returns JSON list of annotation results.
    """
    # ── JSON payload ───────────────────────────────────────────────────────────
    if request.is_json:
        data  = request.get_json(force=True)
        rsids = data.get("rsids", [])

    # ── Form or file upload ────────────────────────────────────────────────────
    else:
        rsid_text = request.form.get("rsid_list", "")
        uploaded  = request.files.get("rsid_file")

        if uploaded and uploaded.filename:
            content   = uploaded.read().decode("utf-8")
            rsid_text = content

        # Split on commas, newlines, tabs, or spaces
        import re
        rsids = [r.strip() for r in re.split(r"[,\n\r\t ]+", rsid_text) if r.strip()]

    if not rsids:
        return jsonify({"error": "No rsIDs provided"}), 400

    if len(rsids) > 200:
        return jsonify({"error": "Maximum 200 rsIDs per batch request"}), 400

    results = batch_analyze_rsids(rsids)

    # Store in session so the download endpoint can re-use without re-computing
    session["last_batch"] = results

    return jsonify({"count": len(results), "results": results})


@app.route("/api/report/<rsid>", methods=["GET"])
def api_report(rsid):
    """Return the full detailed report JSON for a single rsID."""
    report = generate_interpretation_report(rsid.lower())
    return jsonify(report)


# ═══════════════════════════════════════════════════════════════════════════════
# DOWNLOAD ENDPOINT
# ═══════════════════════════════════════════════════════════════════════════════

@app.route("/download/csv", methods=["POST"])
def download_csv():
    """
    Generate and serve a CSV report for the batch analysis results.

    Expects JSON body: { "results": [...] }
    (Re-sends the same results array returned by /api/batch)
    """
    data    = request.get_json(force=True)
    results = data.get("results", [])

    if not results:
        return jsonify({"error": "No results to export"}), 400

    # Build CSV in memory
    output = io.StringIO()
    fieldnames = [
        "rsid", "chromosome", "position", "ref", "alt",
        "gene", "consequence", "hgvs", "protein_change",
        "allele_frequency", "impact_level",
        "clinical_significance", "condition", "pathway",
    ]
    writer = csv.DictWriter(output, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for row in results:
        writer.writerow(row)

    csv_content = output.getvalue()

    response = make_response(csv_content)
    response.headers["Content-Type"] = "text/csv"
    response.headers["Content-Disposition"] = (
        "attachment; filename=variant_analysis_report.csv"
    )
    return response


# ═══════════════════════════════════════════════════════════════════════════════
# ERROR HANDLERS
# ═══════════════════════════════════════════════════════════════════════════════

@app.errorhandler(404)
def not_found(e):
    return render_template("404.html"), 404


@app.errorhandler(500)
def server_error(e):
    return jsonify({"error": "Internal server error", "detail": str(e)}), 500


# ═══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import os

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 5000)))

