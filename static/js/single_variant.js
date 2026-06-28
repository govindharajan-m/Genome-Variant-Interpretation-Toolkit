// ── Mode Tab Switching ────────────────────────────────────────────────────────
const snpTabs = document.querySelectorAll(".mode-tab");
snpTabs.forEach(tab => {
  tab.addEventListener("click", () => {
    snpTabs.forEach(t => { t.classList.remove("active"); t.setAttribute("aria-selected", "false"); });
    tab.classList.add("active");
    tab.setAttribute("aria-selected", "true");
    const mode = tab.dataset.mode;
    document.getElementById("snp-panel-rsid").classList.toggle("hidden", mode !== "rsid");
    document.getElementById("snp-panel-coord").classList.toggle("hidden", mode !== "coord");
  });
});

// ── Coordinate quick-fill ─────────────────────────────────────────────────────
document.querySelectorAll(".example-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.getElementById("chromosome").value = btn.dataset.chr;
    document.getElementById("position").value   = btn.dataset.pos;
    document.getElementById("ref").value        = btn.dataset.ref;
    document.getElementById("alt").value        = btn.dataset.alt;
  });
});

// ── rsID quick-fill ───────────────────────────────────────────────────────────
document.querySelectorAll(".rsid-example-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.getElementById("snp-rsid-input").value = btn.dataset.rsid;
  });
});

// ── rsID Form Submit ──────────────────────────────────────────────────────────
document.getElementById("rsidForm").addEventListener("submit", async e => {
  e.preventDefault();
  const gv-btn  = document.getElementById("rsidBtn");
  const text = document.getElementById("rsidBtnText");
  const spin = document.getElementById("rsidBtnSpinner");

  text.textContent = "Looking up…";
  spin.classList.remove("hidden");
  btn.disabled = true;

  const rsid = document.getElementById("snp-rsid-input").value.trim().toLowerCase();
  if (!rsid) { resetBtn(); return; }

  try {
    const res  = await fetch("/api/analyze-snp-rsid", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({ rsid }),
    });
    
    // We expect the server to return valid JSON even on 404/400
    let data;
    try {
      data = await res.json();
    } catch(e) {
      if (!res.ok) throw new Error("External service temporarily unavailable.");
      throw e;
    }
    
    if (data.error || !data.found) { 
      showSNPError(data.error || "Variant not found."); 
    } else { 
      renderSNPResult(data); 
    }
  } catch(err) {
    showSNPError("Network error — is the server running?");
  } finally {
    text.textContent = "Look Up rsID";
    spin.classList.add("hidden");
    btn.disabled = false;
  }
});

// ── Coordinate Form Submit ────────────────────────────────────────────────────
document.getElementById("snpForm").addEventListener("submit", async e => {
  e.preventDefault();
  const gv-btn  = document.getElementById("analyzeBtn");
  const text = document.getElementById("btnText");
  const spin = document.getElementById("btnSpinner");

  text.textContent = "Analyzing…";
  spin.classList.remove("hidden");
  btn.disabled = true;

  const payload = {
    chromosome: document.getElementById("chromosome").value,
    position:   parseInt(document.getElementById("position").value),
    ref:        document.getElementById("ref").value.toUpperCase(),
    alt:        document.getElementById("alt").value.toUpperCase(),
  };

  try {
    const res  = await fetch("/api/analyze-snp", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify(payload),
    });
    
    let data;
    try {
      data = await res.json();
    } catch(e) {
      if (!res.ok) throw new Error("External service temporarily unavailable.");
      throw e;
    }
    
    if (data.error || !data.found) { 
      showSNPError(data.error || "Variant not found."); 
    } else { 
      renderSNPResult(data); 
    }
  } catch(err) {
    showSNPError("Network error — is the server running?");
  } finally {
    text.textContent = "Analyze Variant";
    spin.classList.add("hidden");
    btn.disabled = false;
  }
});

// ── Render Result ─────────────────────────────────────────────────────────────
function renderSNPResult(d) {
  document.getElementById("emptyState").classList.add("hidden");
  const panel = document.getElementById("resultPanel");
  panel.classList.remove("hidden");

  const isRsid = d.source === "rsid_lookup";

  // ── Header
  const chr = d.input.chromosome;
  const pos = d.input.position;
  document.getElementById("r-variant-title").textContent =
    `Chr${chr}:${pos} ${d.input.ref}>${d.input.alt}`;
  document.getElementById("r-type-badge").textContent = d.classification.type;
  document.getElementById("r-type-badge").className =
    "type-badge type-" + d.classification.type.toLowerCase();

  // rsID badge (top-right of header)
  const rsidBadge = document.getElementById("r-rsid-badge");
  if (isRsid && d.rsid) {
    rsidBadge.textContent = d.rsid;
    rsidBadge.classList.remove("hidden");
  } else {
    rsidBadge.classList.add("hidden");
  }

  // ── Meta row
  document.getElementById("r-chr").textContent     = "Chr " + chr;
  document.getElementById("r-pos").textContent     = pos.toLocaleString();
  document.getElementById("r-alleles").textContent = d.input.ref + " → " + d.input.alt;
  document.getElementById("r-subtype").textContent = d.classification.subtype || "—";

  if (isRsid) {
    document.getElementById("r-hgvs").textContent       = d.hgvs || "—";
    document.getElementById("r-protein").textContent    = d.protein_change || "—";
    const freq = d.allele_frequency;
    document.getElementById("r-allele-freq").textContent =
      freq != null ? (freq * 100).toFixed(2) + "% (" + freq.toFixed(4) + ")" : "—";
  } else {
    document.getElementById("r-hgvs").textContent       = "—";
    document.getElementById("r-protein").textContent    = "—";
    document.getElementById("r-allele-freq").textContent = "—";
  }

  // ── ClinVar Evidence
  if (d.clinvar) {
    document.getElementById("clinvarContainer").innerHTML = renderClinvarHTML(d.clinvar);
  } else {
    document.getElementById("clinvarContainer").innerHTML = "";
  }

  // ── Gene Context & Biological Engine
  if (d.gene_context && d.gene_context.available) {
    document.getElementById("geneContextContainer").innerHTML = renderGeneContextHTML(d.gene_context, d.research_relevance);
    document.getElementById("noGeneCard").classList.add("hidden");
  } else {
    document.getElementById("geneContextContainer").innerHTML = "";
    document.getElementById("noGeneCard").classList.remove("hidden");
  }

  if (d.diseases && d.diseases.available) {
    document.getElementById("diseasesContainer").innerHTML = renderDiseasesHTML(d.diseases);
  } else {
    document.getElementById("diseasesContainer").innerHTML = "";
  }

  if (d.pathways && d.pathways.available) {
    document.getElementById("pathwaysContainer").innerHTML = renderPathwaysHTML(d.pathways);
  } else {
    document.getElementById("pathwaysContainer").innerHTML = "";
  }

  // ── Impact
  const imp      = d.impact;
  const impBadge = document.getElementById("r-impact-badge");
  impBadge.textContent = imp.impact_level;
  impBadge.className   = "impact-badge impact-" + imp.impact_level;
  document.getElementById("r-consequence").textContent =
    imp.consequence.replace(/_/g, " ");
  document.getElementById("r-sift").textContent     = imp.sift_pred || "—";
  document.getElementById("r-polyphen").textContent = imp.polyphen_pred || "—";
  document.getElementById("r-impact-desc").textContent = imp.description;

  // ── Explanation cards (shared renderers from main.js)
  document.getElementById("reportSummaryContainer").innerHTML = renderReportSummaryHTML(d);
  document.getElementById("gwasContainer").innerHTML = renderGWASHTML(d);

  document.getElementById("geneContextContainer").innerHTML = renderGeneContextHTML(d);
  document.getElementById("diseasesContainer").innerHTML = renderDiseasesHTML(d);
  document.getElementById("pathwaysContainer").innerHTML = renderPathwaysHTML(d);

  document.getElementById("populationFrequenciesContainer").innerHTML =
    renderPopulationFrequenciesHTML(d);
  document.getElementById("impactExplanationContainer").innerHTML =
    renderImpactExplanationHTML(d.impact_explanation, imp.impact_level);
  document.getElementById("significanceExplanationContainer").innerHTML =
    renderSignificanceExplanationHTML(
      d.significance_explanation,
      d.clinical_significance || ""
    );
  if (d.interpretation_summary) {
    document.getElementById("evidenceConfidenceContainer").innerHTML = renderEvidenceConfidenceHTML(d);
    document.getElementById("researchRelevanceContainer").innerHTML = renderResearchRelevanceHTML(d);
    document.getElementById("summaryContainer").innerHTML = 
      `<div class="gv-card result-card info-card gv-is-3e33972d" >
         <h3 class="gv-card-title">Summary</h3>
         <p  class="gv-is-9e7fd580">${d.interpretation_summary}</p>
       </div>`;
  } else {
    document.getElementById("evidenceConfidenceContainer").innerHTML = renderEvidenceConfidenceHTML(d);
    document.getElementById("researchRelevanceContainer").innerHTML = renderResearchRelevanceHTML(d);
    document.getElementById("evidenceConfidenceContainer").innerHTML = "";
    document.getElementById("researchRelevanceContainer").innerHTML = "";
    document.getElementById("summaryContainer").innerHTML = "";
  }

  // ── Evidence / references
  document.getElementById("referencesContainer").innerHTML =
    renderEvidenceHTML(d.evidence);
  document.getElementById("pubmedContainer").innerHTML = renderPubMedHTML(d.pubmed);

  
  document.getElementById("dynamicReportNav")?.classList.remove("hidden");
  document.getElementById("dynamicReportTimestamp")?.classList.remove("hidden");

  panel.scrollIntoView({behavior: "smooth"});
}

function showSNPError(msg) {
  document.getElementById("emptyState").classList.add("hidden");
  const panel = document.getElementById("resultPanel");
  panel.classList.remove("hidden");
  
  // Clear all previous results
  document.getElementById("r-variant-title").textContent = "Error";
  document.getElementById("r-type-badge").textContent = "N/A";
  document.getElementById("r-type-badge").className = "type-badge";
  document.getElementById("r-rsid-badge").classList.add("hidden");
  
  document.getElementById("r-chr").textContent = "—";
  document.getElementById("r-pos").textContent = "—";
  document.getElementById("r-alleles").textContent = "—";
  document.getElementById("r-subtype").textContent = "—";
  document.getElementById("r-hgvs").textContent = "—";
  document.getElementById("r-protein").textContent = "—";
  document.getElementById("r-allele-freq").textContent = "—";
  
  document.getElementById("clinvarContainer").innerHTML = "";
  document.getElementById("geneCard").classList.add("hidden");
  document.getElementById("noGeneCard").classList.remove("hidden");
  
  document.getElementById("r-impact-badge").textContent = "—";
  document.getElementById("r-impact-badge").className = "impact-badge";
  document.getElementById("r-consequence").textContent = "—";
  document.getElementById("r-sift").textContent = "—";
  document.getElementById("r-polyphen").textContent = "—";
  document.getElementById("r-impact-desc").textContent = "—";
  
  document.getElementById("gwasContainer").innerHTML = "";
  document.getElementById("populationFrequenciesContainer").innerHTML = "";
  document.getElementById("impactExplanationContainer").innerHTML = "";
  document.getElementById("significanceExplanationContainer").innerHTML = "";
  document.getElementById("evidenceConfidenceContainer").innerHTML = renderEvidenceConfidenceHTML(d);
    document.getElementById("researchRelevanceContainer").innerHTML = renderResearchRelevanceHTML(d);
    document.getElementById("evidenceConfidenceContainer").innerHTML = "";
    document.getElementById("researchRelevanceContainer").innerHTML = "";
    document.getElementById("summaryContainer").innerHTML = "";
  document.getElementById("referencesContainer").innerHTML = "";
  document.getElementById("pubmedContainer").innerHTML = "";

  document.getElementById("dynamicReportNav")?.classList.add("hidden");
  document.getElementById("dynamicReportTimestamp")?.classList.add("hidden");


  // Inject error banner
  const errorHtml = `
    <div class="gv-card result-card info-card gv-is-91ed66af" >
      <h3 class="gv-card-title gv-is-98da24cf" >
        <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
          <circle cx="12" cy="12" r="10"></circle>
          <line x1="12" y1="8" x2="12" y2="12"></line>
          <line x1="12" y1="16" x2="12.01" y2="16"></line>
        </svg>
        Analysis Failed
      </h3>
      <p  class="gv-is-17315441">${msg}</p>
    </div>
  `;
  document.getElementById("referencesContainer").innerHTML = errorHtml;
  
  document.getElementById("dynamicReportNav")?.classList.remove("hidden");
  document.getElementById("dynamicReportTimestamp")?.classList.remove("hidden");

  panel.scrollIntoView({behavior: "smooth"});
}



// Timestamp injector
const tsEl = document.getElementById('dynamicReportTimestampValue');
if (tsEl) tsEl.textContent = new Date().toISOString().replace('T', ' ').substring(0, 16);
