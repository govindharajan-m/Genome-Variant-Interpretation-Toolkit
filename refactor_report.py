import re
with open(r"d:\genome_variant_platform\genome_variant_platform\templates\report.html", "r", encoding="utf-8") as f:
    content = f.read()

# I will extract the blocks from the old content and inject them into a new template.
# The layout should be single-column.

def extract_block(start_marker, end_marker):
    pattern = re.compile(f"{re.escape(start_marker)}(.*?){re.escape(end_marker)}", re.DOTALL)
    match = pattern.search(content)
    if match:
        return match.group(0)
    return ""

variant_overview = extract_block("<!-- Variant Overview -->", "<!-- Functional Impact -->")
functional_impact = extract_block("<!-- Functional Impact -->", "<!-- Impact Explanation -->")
impact_explanation = extract_block("<!-- Impact Explanation -->", "<!-- Clinical Significance Explanation -->")
clinical_significance_explanation = extract_block("<!-- Clinical Significance Explanation -->", "</div><!-- /report-left -->")
clinvar = extract_block("<!-- Clinical Significance -->", "<!-- Population Frequencies -->")
population_frequencies = extract_block("<!-- Population Frequencies -->", "<!-- GWAS Evidence -->")
gwas_evidence = extract_block("<!-- GWAS Evidence -->", "<!-- Version 2.2 Biological Context -->")
gene_context = extract_block("<!-- Version 2.2 Biological Context -->", "</div><!-- /report-right -->")
interpretation_summary = extract_block("<!-- ── Full Interpretation Summary ───────────────────────────────────────── -->", "<!-- ── Evidence & References ─────────────────────────────────────────────── -->")
evidence_references = extract_block("<!-- ── Evidence & References ─────────────────────────────────────────────── -->", "<!-- ── Action Buttons ────────────────────────────────────────────────────── -->")

new_content = content[:content.find('<div class="report-grid">')]

new_content += """
  <!-- 1. Executive Snapshot -->
  <div class="card result-card" style="border-left: 4px solid var(--primary-colour); background-color: var(--bg-secondary);">
    <h2 class="card-title" style="margin-bottom: 12px; font-size: 1.25rem;">EXECUTIVE SNAPSHOT</h2>
    <div class="result-meta-grid" style="row-gap: 1.2rem;">
      <div class="meta-item">
        <span class="meta-label">Variant</span>
        <span class="meta-val mono" style="font-size: 1.1rem; color: var(--text-primary); font-weight: bold;">{{ rsid }}</span>
      </div>
      <div class="meta-item">
        <span class="meta-label">Gene</span>
        <span class="meta-val mono" style="font-size: 1.1rem; color: var(--primary-colour); font-weight: bold;">{% if report.gene_context and report.gene_context.available %}{{ report.gene_context.symbol }}{% else %}Unknown{% endif %}</span>
      </div>
      <div class="meta-item">
        <span class="meta-label">Clinical Significance</span>
        <span class="meta-val" style="font-size: 1.1rem; font-weight: 600;">{{ report.clinvar.clinical_significance if report.clinvar else "Uncertain significance" }}</span>
      </div>
      <div class="meta-item">
        <span class="meta-label">Confidence</span>
        <span class="meta-val" style="font-size: 1.1rem;">{{ report.clinvar.confidence_level if report.clinvar else "Low" }}</span>
      </div>
      <div class="meta-item">
        <span class="meta-label">Research Relevance</span>
        <span class="meta-val" style="font-size: 1.1rem; font-weight: 600; color: #8e44ad;">{{ report.research_relevance or 'Low' }}</span>
      </div>
      <div class="meta-item">
        <span class="meta-label">Population Focus</span>
        <span class="meta-val" style="font-size: 0.95rem; line-height: 1.5;">
          South Asian: <strong>{% if report.population_frequencies and report.population_frequencies.south_asian is not none %}{{ "%.2f"|format(report.population_frequencies.south_asian * 100) }}%{% else %}N/A{% endif %}</strong><br>
          European/Caucasian: <strong>{% if report.population_frequencies and report.population_frequencies.european is not none %}{{ "%.2f"|format(report.population_frequencies.european * 100) }}%{% else %}N/A{% endif %}</strong>
        </span>
      </div>
    </div>
  </div>

  <!-- 2. Interpretation Summary -->
"""

new_content += interpretation_summary + "\n"
new_content += """
  <!-- 3. Variant Overview -->
"""
new_content += variant_overview + "\n"
new_content += functional_impact + "\n"

new_content += """
  <!-- 4. Clinical Evidence -->
"""
new_content += clinvar + "\n"
new_content += clinical_significance_explanation + "\n"

new_content += """
  <!-- 5. Population Analysis -->
"""
new_content += population_frequencies + "\n"

new_content += """
  <!-- 6. Trait & Research Evidence -->
"""
new_content += gwas_evidence + "\n"
new_content += impact_explanation + "\n"

new_content += """
  <!-- 7. Gene Context -->
"""
new_content += gene_context + "\n"

new_content += """
  <!-- 8 & 9. Literature Evidence & References -->
"""
new_content += evidence_references + "\n"

new_content += content[content.find('<!-- ── Action Buttons ────────────────────────────────────────────────────── -->'):]

# Fix up any closing tags that got messed up
new_content = new_content.replace("</div><!-- /report-left -->", "")
new_content = new_content.replace("</div><!-- /report-right -->", "")
new_content = new_content.replace('class="report-right"', 'class="report-content"')
new_content = new_content.replace('class="report-left"', 'class="report-content"')

with open(r"d:\genome_variant_platform\genome_variant_platform\templates\report.html", "w", encoding="utf-8") as f:
    f.write(new_content)
