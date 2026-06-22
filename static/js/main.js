/**
 * main.js — GenomeVAP / Autoradiograph Terminal
 */

function sanitizeCSVValue(value) {
    if (value === null || value === undefined) return '';
    value = String(value);
    if (value.startsWith('=') || value.startsWith('+') || value.startsWith('-') || value.startsWith('@')) {
        return "'" + value;
    }
    return value;
}

document.addEventListener("DOMContentLoaded", () => {

  // ── Auto-uppercase nucleotide inputs ─────────────────────────
  document.querySelectorAll("input[pattern]").forEach(input => {
    input.addEventListener("input", function() {
      const pos = this.selectionStart;
      this.value = this.value.toUpperCase();
      this.setSelectionRange(pos, pos);
    });
  });

  // ── Animate stat numbers on load ─────────────────────────────
  const stats = document.querySelectorAll(".stat-num");
  stats.forEach((el, i) => {
    const original = el.textContent.trim();
    const numMatch = original.match(/[\d.]+/);
    if (!numMatch) return;
    const target = parseFloat(numMatch[0]);
    const prefix = original.slice(0, original.indexOf(numMatch[0]));
    const suffix = original.slice(original.indexOf(numMatch[0]) + numMatch[0].length);
    const isFloat = numMatch[0].includes('.');

    el.textContent = prefix + (isFloat ? '0.0' : '0') + suffix;

    setTimeout(() => {
      const start = performance.now();
      const duration = 900 + i * 120;
      function tick(now) {
        const p = Math.min((now - start) / duration, 1);
        const e = 1 - Math.pow(1 - p, 3);
        const val = isFloat
          ? (target * e).toFixed(1)
          : Math.floor(target * e).toLocaleString();
        el.textContent = prefix + val + suffix;
        if (p < 1) requestAnimationFrame(tick);
        else el.textContent = original;
      }
      requestAnimationFrame(tick);
    }, 200 + i * 80);
  });

  // ── Terminal cursor blink on empty state ─────────────────────
  document.querySelectorAll(".empty-state h3").forEach(el => {
    const orig = el.textContent;
    el.innerHTML = orig + '<span class="cursor">_</span>';
  });
});

/* Utilities */
function formatPos(n) { return n?.toLocaleString() ?? "—"; }
function impactClass(l) { return "impact-" + (l || "MODIFIER"); }

// ── Evidence HTML Generator ──────────────────────────────────────
function renderEvidenceHTML(evidence) {
  if (!evidence) return "";
  
  let html = `<div class="card result-card" id="evidenceCard">
    <h3 class="card-title">Evidence & References</h3>
    <div style="margin-bottom: 1rem; font-weight: bold; color: var(--text-primary);">
      Strength: ${evidence.evidence_strength} (${evidence.evidence_score}/100)
    </div>
    <div class="references-container" style="display: flex; flex-direction: column; gap: 1rem;">`;

  // dbSNP
  html += `<div class="ref-section">
    <div style="font-weight: 600; margin-bottom: 0.25rem; color: var(--text-primary);">dbSNP</div>`;
  if (evidence.dbsnp && evidence.dbsnp.id) {
    html += `<a href="${evidence.dbsnp.url}" target="_blank" class="ext-link" style="color: var(--primary-colour); text-decoration: none;">${evidence.dbsnp.id}</a>`;
  } else {
    html += `<div style="color: var(--text-secondary); font-size: 0.875rem;">No reference available</div>`;
  }
  html += `</div>`;

  // ClinVar
  html += `<div class="ref-section">
    <div style="font-weight: 600; margin-bottom: 0.25rem; color: var(--text-primary);">ClinVar</div>`;
  if (evidence.clinvar && evidence.clinvar.id) {
    html += `<div>
      <a href="${evidence.clinvar.url}" target="_blank" class="ext-link" style="color: var(--primary-colour); text-decoration: none;">${evidence.clinvar.id}</a>
      <div style="font-size: 0.875rem; color: var(--text-secondary); margin-top: 0.25rem;">
        Classification: ${evidence.clinvar.classification} <br/>
        Review Status: ${evidence.clinvar.review_status}
      </div>
    </div>`;
  } else {
    html += `<div style="color: var(--text-secondary); font-size: 0.875rem;">No reference available</div>`;
  }
  html += `</div>`;

  // GeneCards
  html += `<div class="ref-section">
    <div style="font-weight: 600; margin-bottom: 0.25rem; color: var(--text-primary);">GeneCards</div>`;
  if (evidence.genecards && evidence.genecards.id) {
    html += `<a href="${evidence.genecards.url}" target="_blank" class="ext-link" style="color: var(--primary-colour); text-decoration: none;">${evidence.genecards.id}</a>`;
  } else {
    html += `<div style="color: var(--text-secondary); font-size: 0.875rem;">No reference available</div>`;
  }
  html += `</div>`;

  // PubMed
  html += `<div class="ref-section">
    <div style="font-weight: 600; margin-bottom: 0.25rem; color: var(--text-primary);">PubMed</div>`;
  if (evidence.pubmed && evidence.pubmed.length > 0) {
    html += `<ul style="list-style-type: disc; margin: 0; padding-left: 1.5rem;">`;
    for (const item of evidence.pubmed) {
      html += `<li style="margin-bottom: 0.5rem; font-size: 0.875rem;">
        <a href="${item.url}" target="_blank" class="ext-link" style="color: var(--primary-colour); text-decoration: none;">
          ${item.id}
        </a>: ${item.title} (${item.year})
      </li>`;
    }
    html += `</ul>`;
  } else {
    html += `<div style="color: var(--text-secondary); font-size: 0.875rem;">No reference available</div>`;
  }
  html += `</div>`;

  if (evidence.warnings && evidence.warnings.length > 0) {
    html += `<div style="margin-top: 1rem; color: #e74c3c; font-size: 0.875rem;">
      <strong>Warnings:</strong>
      <ul style="margin: 0; padding-left: 1.5rem;">`;
    for (const w of evidence.warnings) {
      html += `<li>${w}</li>`;
    }
    html += `</ul></div>`;
  }

  html += `</div></div>`;
  return html;
}

// ── Impact Explanation HTML Generator ────────────────────────────────────────
function renderImpactExplanationHTML(impactExplanation, impactLevel) {
  if (!impactExplanation) return "";

  const level = (impactLevel || "MODIFIER").toUpperCase();
  const impactColours = {
    HIGH: "var(--high, #c84848)",
    MODERATE: "var(--moderate, #c07830)",
    LOW: "var(--low, #4a9868)",
    MODIFIER: "var(--modifier, #4e6878)"
  };
  const colour = impactColours[level] || impactColours.MODIFIER;

  let consequencesHtml = "";
  if (impactExplanation.biological_consequences && impactExplanation.biological_consequences.length) {
    consequencesHtml = `<ul class="impact-expl-list">` +
      impactExplanation.biological_consequences.map(c => `<li>${c}</li>`).join("") +
      `</ul>`;
  }

  let examplesHtml = "";
  if (impactExplanation.example_variant_types && impactExplanation.example_variant_types.length) {
    examplesHtml = impactExplanation.example_variant_types
      .map(e => `<span class="impact-expl-tag">${e}</span>`)
      .join("");
  }

  return `<div class="card result-card impact-expl-card" id="impactExplanationCard">
    <details class="impact-expl-details">
      <summary class="impact-expl-summary">
        <span class="impact-expl-summary-text">
          <span class="impact-expl-icon">?</span>
          Impact Explanation — What does <span class="impact-expl-level" style="color: ${colour};">${level}</span> mean?
        </span>
        <span class="impact-expl-toggle">▸</span>
      </summary>
      <div class="impact-expl-body">
        <div class="impact-expl-section">
          <div class="impact-expl-label">Meaning</div>
          <p class="impact-expl-text">${impactExplanation.meaning || "—"}</p>
        </div>
        <div class="impact-expl-section">
          <div class="impact-expl-label">Typical Biological Consequences</div>
          ${consequencesHtml}
        </div>
        <div class="impact-expl-section">
          <div class="impact-expl-label">Example Variant Types</div>
          <div class="impact-expl-tags">${examplesHtml}</div>
        </div>
        <div class="impact-expl-section">
          <div class="impact-expl-label">Interpretation Guidance</div>
          <p class="impact-expl-interp">${impactExplanation.interpretation || "—"}</p>
        </div>
      </div>
    </details>
  </div>`;
}

// ── Clinical Significance Explanation HTML Generator ─────────────────────────
function renderSignificanceExplanationHTML(sigExplanation, significance) {
  if (!sigExplanation) return "";

  const sig = (significance || "Unknown").replace(/_/g, " ");
  const sigColours = {
    "Pathogenic":            "#e74c3c",
    "Likely_pathogenic":     "#e67e22",
    "Uncertain_significance":"#f39c12",
    "Likely_benign":         "#2ecc71",
    "Benign":                "#27ae60",
    "Risk_factor":           "#9b59b6",
    "drug_response":         "#3498db",
  };
  const colour = sigColours[significance] || "var(--amber)";

  let evidenceHtml = "";
  if (sigExplanation.evidence_basis && sigExplanation.evidence_basis.length) {
    evidenceHtml = `<ul class="sig-expl-list">` +
      sigExplanation.evidence_basis.map(e => `<li>${e}</li>`).join("") +
      `</ul>`;
  }

  return `<div class="card result-card sig-expl-card" id="significanceExplanationCard">
    <details class="sig-expl-details">
      <summary class="sig-expl-summary">
        <span class="sig-expl-summary-text">
          <span class="sig-expl-icon">⚕</span>
          Clinical Significance Explanation — <span class="sig-expl-level" style="color: ${colour};">${sig}</span>
        </span>
        <span class="sig-expl-toggle">▸</span>
      </summary>
      <div class="sig-expl-body">
        <div class="sig-expl-section">
          <div class="sig-expl-label">Definition</div>
          <p class="sig-expl-text">${sigExplanation.definition || "—"}</p>
        </div>
        <div class="sig-expl-section">
          <div class="sig-expl-label">Interpretation</div>
          <p class="sig-expl-interp">${sigExplanation.interpretation || "—"}</p>
        </div>
        <div class="sig-expl-section">
          <div class="sig-expl-label">Evidence Basis</div>
          ${evidenceHtml}
        </div>
      </div>
    </details>
  </div>`;
}

// ── Population Frequencies HTML Generator ───────────────────────────────────
function renderPopulationFrequenciesHTML(data) {
  if (!data || !data.population_frequencies) return "";

  const freqs = data.population_frequencies;
  const interp = data.frequency_interpretation || "";
  const source = data.frequency_source || "";
  const comp = data.frequency_comparison;

  const renderRow = (label, key, isGlobal = false) => {
    const val = freqs[key];
    const width = val != null ? Math.min(val * 100, 100) : 0;
    const textVal = val != null ? (val * 100).toFixed(2) + "%" : "—";
    const barClass = isGlobal ? "pop-freq-bar global" : "pop-freq-bar";
    return `
      <div class="pop-freq-row">
        <span class="pop-freq-label">${label}</span>
        <div class="pop-freq-bar-container"><div class="${barClass}" style="width: ${width}%;"></div></div>
        <span class="pop-freq-val">${textVal}</span>
      </div>
    `;
  };

  let sourceHtml = "";
  if (source) {
    sourceHtml = `<div class="pop-freq-source">Source: ${source}</div>`;
  }
  
  let comparisonHtml = "";
  if (comp) {
    let foldText = comp.fold_difference === Infinity ? "Infinity" : comp.fold_difference.toFixed(1) + "×";
    comparisonHtml = `
      <div style="margin-top: 1rem; padding: 1rem; background: var(--bg-secondary); border-radius: 0.5rem; border-left: 3px solid var(--primary-colour);">
        <h4 style="margin-top: 0; margin-bottom: 0.5rem; font-size: 0.95rem; color: var(--text-primary);">Population Comparison</h4>
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem;">
          <div>
            <span style="font-size: 0.85rem; color: var(--text-secondary); display: block;">Absolute Difference</span>
            <span style="font-weight: 600; font-size: 1.1rem; color: var(--text-primary);">${(comp.absolute_difference * 100).toFixed(2)}%</span>
          </div>
          <div>
            <span style="font-size: 0.85rem; color: var(--text-secondary); display: block;">Fold Difference</span>
            <span style="font-weight: 600; font-size: 1.1rem; color: var(--text-primary);">${foldText}</span>
          </div>
        </div>
      </div>
    `;
  }

  return `
  <div class="card result-card" id="populationFrequenciesCard">
    <h3 class="card-title">Population Frequencies</h3>
    <p class="interpretation-text" style="margin-top: 0;">${interp}</p>
    <div class="pop-freq-grid">
      ${renderRow("Global", "global", true)}
      ${renderRow("South Asian", "south_asian")}
      ${renderRow("European", "european")}
    </div>
    ${comparisonHtml}
    ${sourceHtml}
  </div>`;
}

// ── GWAS Evidence HTML Generator ───────────────────────────────────────────
function renderGWASHTML(data) {
  if (!data) return "";
  
  const topHit = data.gwas_top_hit;
  const interp = data.gwas_interpretation;

  if (!topHit) {
    return `
    <div class="card result-card" id="gwasCard">
      <h3 class="card-title">GWAS Evidence</h3>
      <div class="gwas-empty">No GWAS associations available.</div>
    </div>`;
  }

  const formatSci = (val) => {
    if (val == null) return "—";
    if (val < 0.001 || val > 1000) return val.toExponential(2);
    return val.toString();
  };
  
  const studyLink = topHit.pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${topHit.pmid}/" target="_blank">${topHit.study || topHit.pmid} ↗</a>` : (topHit.study || "—");

  return `
  <div class="card result-card" id="gwasCard">
    <h3 class="card-title">GWAS Evidence</h3>
    <div class="gwas-meta-grid">
      <div class="gwas-meta-item">
        <span class="gwas-meta-label">Trait</span>
        <span class="gwas-meta-val">${topHit.trait || "—"}</span>
      </div>
      <div class="gwas-meta-item">
        <span class="gwas-meta-label">P-value</span>
        <span class="gwas-meta-val mono">${formatSci(topHit.p_value)}</span>
      </div>
      <div class="gwas-meta-item">
        <span class="gwas-meta-label">Odds Ratio</span>
        <span class="gwas-meta-val mono">${topHit.odds_ratio != null ? topHit.odds_ratio.toFixed(2) : "—"}</span>
      </div>
      <div class="gwas-meta-item">
        <span class="gwas-meta-label">Beta</span>
        <span class="gwas-meta-val mono">${topHit.beta != null ? topHit.beta.toFixed(3) : "—"}</span>
      </div>
      <div class="gwas-meta-item">
        <span class="gwas-meta-label">Study</span>
        <span class="gwas-meta-val">${studyLink}</span>
      </div>
    </div>
    
    ${interp ? `
    <div class="gwas-interp">
      <strong>${interp.association_strength}</strong>: ${interp.narrative}
    </div>
    ` : ""}
    
    <div class="gwas-source">Source: GWAS Catalog REST API</div>
  </div>`;
}
// ── Report Summary Generator ───────────────────────────────────────────
function renderReportSummaryHTML(d) {
  if (!d) return "";
  
  // Variant
  const variant = d.rsid || (d.input ? `Chr${d.input.chromosome}:${d.input.position || d.input.start}` : "Unknown");
  
  // Gene
  let gene = "Unknown";
  if (d.gene_context && d.gene_context.available) {
      gene = d.gene_context.symbol;
  } else if (d.gene && d.gene.gene) {
      gene = d.gene.gene;
  }
  
  // Clinical Significance & Confidence
  let clinSig = "Uncertain significance";
  let clinConf = "Low";
  if (d.clinvar) {
      clinSig = d.clinvar.clinical_significance || clinSig;
      clinConf = d.clinvar.confidence_level || clinConf;
  }
  
  // Research Relevance
  const relevance = (typeof d.research_relevance === "object") ? (d.research_relevance.tier || "Low") : (d.research_relevance || "Low");
  
  // Populations
  let sas = "N/A";
  let eur = "N/A";
  if (d.population_frequencies) {
      if (d.population_frequencies.South_Asian != null) {
          sas = (d.population_frequencies.South_Asian * 100).toFixed(2) + "%";
      }
      if (d.population_frequencies.European != null) {
          eur = (d.population_frequencies.European * 100).toFixed(2) + "%";
      }
  }

  return `
    <div class="card result-card" style="border-left: 4px solid var(--primary-colour); background-color: var(--bg-secondary);">
      <h2 class="card-title" style="margin-bottom: 12px; font-size: 1.25rem;">REPORT SUMMARY</h2>
      <div class="result-meta-grid" style="row-gap: 1.2rem;">
        <div class="meta-item">
          <span class="meta-label">Variant</span>
          <span class="meta-val mono" style="font-size: 1.1rem; color: var(--text-primary); font-weight: bold;">${variant}</span>
        </div>
        <div class="meta-item">
          <span class="meta-label">Gene</span>
          <span class="meta-val mono" style="font-size: 1.1rem; color: var(--primary-colour); font-weight: bold;">${gene}</span>
        </div>
        <div class="meta-item">
          <span class="meta-label">Clinical Significance</span>
          <span class="meta-val" style="font-size: 1.1rem; font-weight: 600;">${clinSig}</span>
        </div>
        <div class="meta-item">
          <span class="meta-label">Confidence</span>
          <span class="meta-val" style="font-size: 1.1rem;">${clinConf}</span>
        </div>
        <div class="meta-item">
          <span class="meta-label">Research Relevance</span>
          <span class="meta-val" style="font-size: 1.1rem; font-weight: 600;">${relevance}</span>
        </div>
        <div class="meta-item">
          <span class="meta-label">Population Focus</span>
          <span class="meta-val" style="font-size: 0.95rem; line-height: 1.5;">
            South Asian: <strong>${sas}</strong><br>
            European/Caucasian: <strong>${eur}</strong>
          </span>
        </div>
      </div>
    </div>
  `;
}

// ── Version 2.2 Biological Context Generators ───────────────────────────────
function renderGeneContextHTML(data, relevance) {
  if (!data || !data.available) return "";
  
  let relColor = "#7f8c8d";
  if (relevance === "Very High") relColor = "#8e44ad";
  else if (relevance === "High") relColor = "#2980b9";
  else if (relevance === "Moderate") relColor = "#27ae60";

  return `
    <div class="card result-card" style="margin-top: 1.2rem; padding-top: 1rem; border-top: 1px solid var(--border);">
      <h3 class="card-title">Gene Context</h3>
      <div class="result-meta-grid" style="margin-bottom: 1rem;">
        <div class="meta-item"><span class="meta-label">Gene Symbol</span><span class="meta-val mono" style="font-weight:bold; color:var(--primary-colour);">${data.symbol}</span></div>
        <div class="meta-item"><span class="meta-label">Full Name</span><span class="meta-val">${data.full_name || "—"}</span></div>
        <div class="meta-item"><span class="meta-label">Chromosome</span><span class="meta-val mono">${data.location || "—"}</span></div>
        <div class="meta-item">
          <span class="meta-label">Research Relevance</span>
          <span class="meta-val" style="font-weight:600; color:${relColor};">${relevance || "Low"}</span>
        </div>
      </div>
      <div class="gene-description">
        <span style="font-weight: 600; font-size: 0.9rem; color: var(--text-secondary); display: block; margin-bottom: 0.25rem;">Function:</span>
        ${data.description || "No official description available."}
      </div>
    </div>
  `;
}

function renderDiseasesHTML(data) {
  if (!data || !data.available || !data.diseases.length) return "";
  const listItems = data.diseases.map(d => `<li>${d}</li>`).join("");
  return `
    <div class="card result-card">
      <h3 class="card-title">Associated Diseases</h3>
      <ul style="margin: 0; padding-left: 1.5rem; color: var(--text-primary); line-height: 1.6;">
        ${listItems}
      </ul>
    </div>
  `;
}

function renderPathwaysHTML(data) {
  if (!data || !data.available || !data.pathways.length) return "";
  const listItems = data.pathways.map(p => `<li>${p}</li>`).join("");
  return `
    <div class="card result-card">
      <h3 class="card-title">Biological Pathways</h3>
      <ul style="margin: 0; padding-left: 1.5rem; color: var(--text-primary); line-height: 1.6;">
        ${listItems}
      </ul>
    </div>
  `;
}


// ── ClinVar Evidence HTML Generator ─────────────────────────────────────────
function renderClinvarHTML(data) {
  if (!data) return `
    <div class="card result-card" id="clinvarCard">
      <h3 class="card-title">ClinVar Evidence</h3>
      <div class="gwas-empty">No ClinVar annotation available.</div>
    </div>`;

  const confClassMap = {
    "Very High": "cv-conf-very-high",
    "High": "cv-conf-high",
    "Moderate": "cv-conf-moderate",
    "Low": "cv-conf-low"
  };
  
  const confClass = confClassMap[data.confidence_level] || "cv-conf-low";
  
  const sigColors = {
    "Pathogenic":"#e74c3c","Likely pathogenic":"#e67e22","Likely_pathogenic":"#e67e22",
    "Uncertain significance":"#f39c12","Uncertain_significance":"#f39c12",
    "Likely benign":"#2ecc71","Likely_benign":"#2ecc71",
    "Benign":"#27ae60","Risk factor":"#9b59b6","Drug response":"#3498db"
  };
  
  // Normalize key
  let sigKey = data.clinical_significance;
  if (sigKey) sigKey = sigKey.charAt(0).toUpperCase() + sigKey.slice(1);
  const color = sigColors[sigKey] || "var(--text-primary)";

  let conflictHtml = "";
  if (data.conflicting_interpretations && data.conflict_breakdown && Object.keys(data.conflict_breakdown).length > 0) {
    let breakdownHtml = "";
    for (const [cls, count] of Object.entries(data.conflict_breakdown)) {
      breakdownHtml += `
        <div class="cv-conflict-item">
          <span>${cls}</span>
          <span class="cv-count">${count}</span>
        </div>`;
    }
    
    conflictHtml = `
      <div class="cv-conflict-panel">
        <div class="cv-conflict-title">
          Conflicting Interpretations
          <span class="cv-consensus">Consensus: ${data.consensus_score != null ? data.consensus_score + '%' : 'N/A'}</span>
        </div>
        <div class="cv-conflict-list">
          ${breakdownHtml}
        </div>
      </div>
    `;
  }

  return `
  <div class="card result-card" id="clinvarCard">
    <div class="result-header">
      <div>
        <h3 class="card-title" style="margin-bottom: 0;">ClinVar Evidence</h3>
      </div>
      <div>
        <span class="status-badge badge-high">${data.confidence_level} Confidence</span>
      </div>
    </div>
    
    <div class="result-meta-grid" style="margin-top: 1rem;">
      <div class="meta-item">
        <span class="meta-label">Clinical Significance</span>
        <span class="status-badge badge-pathogenic" style="color: ${color}; border-color: ${color};">${data.clinical_significance || "—"}</span>
      </div>
      <div class="meta-item">
        <span class="meta-label">Condition</span>
        <span class="meta-val">${data.condition || "—"}</span>
      </div>
      <div class="meta-item">
        <span class="meta-label">Review Status</span>
        <span class="meta-val">${data.review_status || "—"}</span>
      </div>
      <div class="meta-item">
        <span class="meta-label">Accession</span>
        <span class="meta-val mono">
          <a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${data.variation_id}/" target="_blank" style="color: var(--accent); text-decoration: none;">
            ${data.accession || data.variation_id} ↗
          </a>
        </span>
      </div>
      <div class="meta-item">
        <span class="meta-label">Last Evaluated</span>
        <span class="meta-val mono">${data.last_evaluated || "—"}</span>
      </div>
    </div>
    
    ${conflictHtml}
  </div>`;
}

function renderPubMedHTML(data) {
  if (!data || !data.available || data.paper_count === 0) {
    return `
    <div class="card result-card" style="margin-top: 1rem;">
      <div class="result-header">
        <h3 class="card-title" style="margin-bottom: 0;">Literature Evidence</h3>
      </div>
      <p style="margin-top: 1rem; color: var(--text-secondary);">No PubMed evidence available.</p>
    </div>
    `;
  }

  let papersHtml = "";
  data.papers.forEach((p, idx) => {
    papersHtml += `
      <div style="padding: 0.75rem 0; border-bottom: 1px solid var(--border-2); ${idx === data.papers.length - 1 ? 'border-bottom: none;' : ''}">
        <div style="font-weight: 500; margin-bottom: 0.25rem;">
          <a href="${p.url}" target="_blank" style="color: var(--accent); text-decoration: none;">${p.title}</a>
        </div>
        <div style="font-size: 0.85rem; color: var(--text-secondary);">
          ${p.journal} (${p.year}) &bull; PMID: <a href="${p.url}" target="_blank" style="color: var(--text-secondary);">${p.pmid}</a> &bull; ${p.has_abstract ? "Abstract available" : "No abstract"}
        </div>
      </div>
    `;
  });

  return `
  <div class="card result-card" style="margin-top: 1rem;">
    <div class="result-header" style="align-items: flex-start;">
      <div>
        <h3 class="card-title" style="margin-bottom: 0;">Literature Evidence</h3>
        <div style="font-size: 0.9rem; color: var(--text-secondary); margin-top: 0.25rem;">
          ${data.paper_count} Publication${data.paper_count > 1 ? 's' : ''} Found
        </div>
      </div>
      <div>
        <span class="type-badge" style="background: rgba(43, 144, 217, 0.1); color: var(--accent); border: 1px solid rgba(43, 144, 217, 0.2); font-size: 0.85rem; padding: 0.3rem 0.6rem;">
          ${data.evidence_level || "Unknown Evidence"}
        </span>
      </div>
    </div>
    
    <div style="margin-top: 1rem;">
      <h4 style="font-size: 0.95rem; text-transform: uppercase; letter-spacing: 0.5px; color: var(--text-secondary); margin-bottom: 0.5rem; border-bottom: 1px solid var(--border-2); padding-bottom: 0.25rem;">Top Publications</h4>
      ${papersHtml}
    </div>
  </div>
  `;
}


// ── Evidence Confidence HTML Generator ─────────────────────────────────────
function renderEvidenceConfidenceHTML(d) {
    if (!d || !d.evidence_confidence) return "";
    const ev = d.evidence_confidence;
    
    let factorsHtml = "";
    if (ev.contributing_factors && ev.contributing_factors.length > 0) {
        factorsHtml = `<div style="margin-top: 1rem;"><h4 style="font-size: 0.95rem; margin-bottom: 0.5rem; color: var(--text-primary);">Contributing Factors</h4><ul style="list-style-type: none; padding: 0; margin: 0;">`;
        ev.contributing_factors.forEach(f => {
            factorsHtml += `<li style="font-size: 0.9rem; color: var(--text-secondary); margin-bottom: 0.3rem;">✓ ${f}</li>`;
        });
        factorsHtml += `</ul></div>`;
    }
    
    let limitsHtml = "";
    if (ev.limitations && ev.limitations.length > 0) {
        limitsHtml = `<div style="margin-top: 1rem;"><h4 style="font-size: 0.95rem; margin-bottom: 0.5rem; color: var(--text-primary);">Limitations</h4><ul style="list-style-type: none; padding: 0; margin: 0;">`;
        ev.limitations.forEach(l => {
            limitsHtml += `<li style="font-size: 0.9rem; color: var(--text-secondary); margin-bottom: 0.3rem;">• ${l}</li>`;
        });
        limitsHtml += `</ul></div>`;
    }

    return `
    <div class="card result-card" id="evidenceConfidenceCard" style="margin-top: 1rem;">
        <h3 class="card-title">Evidence Confidence</h3>
        <div style="display: flex; align-items: baseline; gap: 1rem; margin-bottom: 1rem;">
            <span style="font-size: 2rem; font-weight: bold; color: var(--primary-colour);">${ev.score} <span style="font-size: 1rem; color: var(--text-secondary);">/ 100</span></span>
            <span class="status-badge badge-strong" style="font-size: 1rem; padding: 0.4rem 0.8rem; background: var(--bg-tertiary);">${ev.tier.toUpperCase()} CONFIDENCE</span>
            <span style="font-weight: 600; color: var(--text-primary);">${ev.strength}</span>
        </div>
        <p style="font-size: 0.95rem; line-height: 1.6; color: var(--text-secondary);">${ev.narrative}</p>
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem;">
            ${factorsHtml}
            ${limitsHtml}
        </div>
    </div>`;
}


// ── Research Relevance HTML Generator ─────────────────────────────────────
function renderResearchRelevanceHTML(d) {
    if (!d || !d.research_relevance || typeof d.research_relevance !== 'object') return "";
    const rr = d.research_relevance;
    
    let reasonsHtml = "";
    if (rr.reasons && rr.reasons.length > 0) {
        reasonsHtml = `<div style="margin-top: 1rem;"><h4 style="font-size: 0.95rem; margin-bottom: 0.5rem; color: var(--text-primary);">Why This Variant Matters</h4><ul style="list-style-type: none; padding: 0; margin: 0;">`;
        rr.reasons.forEach(r => {
            reasonsHtml += `<li style="font-size: 0.9rem; color: var(--text-secondary); margin-bottom: 0.3rem;">✓ ${r}</li>`;
        });
        reasonsHtml += `</ul></div>`;
    }
    
    let appsHtml = "";
    if (rr.applications && rr.applications.length > 0) {
        appsHtml = `<div style="margin-top: 1rem;"><h4 style="font-size: 0.95rem; margin-bottom: 0.5rem; color: var(--text-primary);">Potential Applications</h4><ul style="list-style-type: none; padding: 0; margin: 0;">`;
        rr.applications.forEach(a => {
            appsHtml += `<li style="font-size: 0.9rem; color: var(--text-secondary); margin-bottom: 0.3rem;">• ${a}</li>`;
        });
        appsHtml += `</ul></div>`;
    }

    return `
    <div class="card result-card" id="researchRelevanceCard" style="margin-top: 1rem;">
        <h3 class="card-title">Research Relevance</h3>
        <div style="display: flex; align-items: baseline; gap: 1rem; margin-bottom: 1rem;">
            <span class="status-badge badge-strong" style="font-size: 1rem; padding: 0.4rem 0.8rem; background: var(--bg-tertiary);">${rr.tier ? rr.tier.toUpperCase() : "UNKNOWN"}</span>
            <span style="font-size: 1rem; color: var(--text-secondary);">Research Score: <span style="font-weight: 600; color: var(--primary-colour);">${rr.score} / 100</span></span>
        </div>
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem;">
            ${reasonsHtml}
            ${appsHtml}
        </div>
    </div>`;
}


// ── ACMG Evidence HTML Generator ─────────────────────────────────────
function renderACMGEvidenceHTML(d) {
    if (!d || !d.acmg_evidence) return "";
    const acmg = d.acmg_evidence;
    
    let criteriaHtml = "";
    if (acmg.criteria && acmg.criteria.length > 0) {
        criteriaHtml = `<div style="margin-top: 1rem;"><h4 style="font-size: 0.95rem; margin-bottom: 0.5rem; color: var(--text-primary);">Mapped Criteria</h4><ul style="list-style-type: none; padding: 0; margin: 0;">`;
        acmg.criteria.forEach(c => {
            const codeColor = c.code.startsWith('P') ? '#e74c3c' : '#27ae60';
            criteriaHtml += `<li style="font-size: 0.95rem; color: var(--text-secondary); margin-bottom: 0.8rem;">
                <span style="font-weight: 600; color: ${codeColor}; border: 1px solid ${codeColor}; padding: 0.1rem 0.3rem; border-radius: 3px; font-size: 0.8rem;">${c.code}</span>
                <span style="font-weight: 600; margin-left: 0.3rem;">— ${c.strength}</span><br>
                <span style="margin-left: 2rem; display: block; font-size: 0.9rem; margin-top: 0.2rem;">${c.reason}</span>
            </li>`;
        });
        criteriaHtml += `</ul></div>`;
    }

    return `
    <div class="card result-card" id="acmgEvidenceCard" style="margin-top: 1rem; border-left: 4px solid #9b59b6;">
        <h3 class="card-title">ACMG Evidence Mapping</h3>
        
        <div style="background-color: var(--bg-secondary); padding: 1rem; border-radius: 6px; margin-bottom: 1.5rem; border-left: 3px solid #f39c12;">
            <p style="margin: 0; font-size: 0.9rem; color: var(--text-secondary);">
                <strong>Research Use Only:</strong> This module maps potentially relevant ACMG evidence categories based on available platform data. It is not a clinical ACMG classifier and must not be used for diagnosis, treatment decisions, or patient management.
            </p>
        </div>

        <div style="display: flex; align-items: baseline; gap: 1rem; margin-bottom: 1rem;">
            <span class="status-badge badge-strong" style="font-size: 1.1rem; padding: 0.4rem 0.8rem; background: var(--bg-tertiary);">${acmg.classification.toUpperCase()}</span>
            <span style="font-size: 0.9rem; color: var(--text-secondary);">Mapping Confidence: <strong style="color: var(--text-primary);">${acmg.mapping_confidence}</strong></span>
        </div>
        
        <div style="display: flex; gap: 2rem; margin-bottom: 1rem;">
            <span style="font-size: 0.9rem; color: var(--text-secondary);">Pathogenic Triggers: <strong style="color: var(--text-primary);">${acmg.pathogenic_evidence_count}</strong></span>
            <span style="font-size: 0.9rem; color: var(--text-secondary);">Benign Triggers: <strong style="color: var(--text-primary);">${acmg.benign_evidence_count}</strong></span>
        </div>
        
        ${criteriaHtml}
    </div>`;
}


// ── Pharmacogenomics HTML Generator ─────────────────────────────────────
function renderPharmacogenomicsHTML(d) {
    if (!d || !d.pharmacogenomics || d.pharmacogenomics.tier === "Not Available") return "";
    const pgx = d.pharmacogenomics;
    
    let appsHtml = "";
    if (pgx.applications && pgx.applications.length > 0) {
        appsHtml = `<div style="margin-top: 1rem;"><h4 style="font-size: 0.95rem; margin-bottom: 0.5rem; color: var(--text-primary);">Applications</h4><ul style="list-style-type: none; padding: 0; margin: 0;">`;
        pgx.applications.forEach(a => {
            appsHtml += `<li style="font-size: 0.9rem; color: var(--text-secondary); margin-bottom: 0.3rem;">✓ ${a}</li>`;
        });
        appsHtml += `</ul></div>`;
    }

    let interactionsHtml = "";
    if (pgx.interactions && pgx.interactions.length > 0) {
        interactionsHtml = `<div style="margin-top: 1rem;"><h4 style="font-size: 0.95rem; margin-bottom: 0.5rem; color: var(--text-primary);">Drug Relevance</h4><ul style="list-style-type: none; padding: 0; margin: 0;">`;
        pgx.interactions.forEach(i => {
            interactionsHtml += `<li style="font-size: 0.95rem; color: var(--text-secondary); margin-bottom: 0.8rem; background: var(--bg-secondary); padding: 0.8rem; border-radius: 4px; border-left: 3px solid #16a085;">
                <span style="font-weight: bold; color: var(--primary-colour); font-size: 1.05rem;">${i.drug}</span>
                <span class="status-badge badge-strong" style="font-size: 0.75rem; margin-left: 0.5rem;">${i.evidence_level} Evidence</span><br>
                <span style="display: block; margin-top: 0.4rem; font-size: 0.9rem;">${i.association}</span>
            </li>`;
        });
        interactionsHtml += `</ul></div>`;
    }

    return `
    <div class="card result-card" id="pharmacogenomicsCard" style="margin-top: 1rem; border-left: 4px solid #16a085;">
        <h3 class="card-title">Pharmacogenomics & Precision Medicine</h3>
        
        <div style="background-color: var(--bg-secondary); padding: 1rem; border-radius: 6px; margin-bottom: 1.5rem; border-left: 3px solid #f39c12;">
            <p style="margin: 0; font-size: 0.9rem; color: var(--text-secondary);">
                <strong>Research Use Only:</strong> This module identifies potential pharmacogenomic and precision medicine relevance based on curated databases (e.g., PharmGKB, CPIC). It must not be used to recommend medications, suggest doses, or provide treatment advice.
            </p>
        </div>

        <div style="display: flex; align-items: baseline; gap: 1rem; margin-bottom: 1rem;">
            <span class="status-badge badge-strong" style="font-size: 1.1rem; padding: 0.4rem 0.8rem; background: var(--bg-tertiary);">${pgx.tier.toUpperCase()} TIER</span>
            <span style="font-size: 0.9rem; color: var(--text-secondary);">Relevance Score: <strong style="color: var(--text-primary);">${pgx.score}</strong>/100</span>
        </div>
        
        ${interactionsHtml}
        ${appsHtml}
        
        <div style="margin-top: 1.5rem; font-size: 0.85rem; color: var(--text-secondary);">
            <strong>Evidence Source:</strong> Curated mappings derived from authoritative pharmacogenomic resources (PharmGKB / CPIC / FDA).
        </div>
    </div>`;
}
