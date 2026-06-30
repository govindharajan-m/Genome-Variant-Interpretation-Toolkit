document.addEventListener("DOMContentLoaded", () => {
// ── Mode Tab Switching ────────────────────────────────────────────────────────
const cnvTabs = document.querySelectorAll(".mode-tab");
cnvTabs.forEach(tab => {
  tab.addEventListener("click", () => {
    cnvTabs.forEach(t => { t.classList.remove("active"); t.setAttribute("aria-selected", "false"); });
    tab.classList.add("active");
    tab.setAttribute("aria-selected", "true");
    const mode = tab.dataset.mode;
    document.getElementById("cnv-panel-rsid").classList.toggle("hidden", mode !== "rsid");
    document.getElementById("cnv-panel-gene").classList.toggle("hidden", mode !== "gene");
    document.getElementById("cnv-panel-coord").classList.toggle("hidden", mode !== "coord");
  });
});

// ── CNV Type Toggle ───────────────────────────────────────────────────────────
document.querySelectorAll(".cnv-toggle").forEach(btn => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".cnv-toggle").forEach(b => b.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById("cnv-type").value = btn.dataset.type;
  });
});

// ── rsID quick-fill ───────────────────────────────────────────────────────────
document.querySelectorAll(".cnv-rsid-example-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.getElementById("cnv-rsid-input").value = btn.dataset.rsid;
  });
});

// ── Gene quick-fill ───────────────────────────────────────────────────────────
document.querySelectorAll(".cnv-gene-example-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.getElementById("cnv-gene-input").value = btn.dataset.gene;
  });
});

// ── Coordinate quick-fill ─────────────────────────────────────────────────────
document.querySelectorAll(".cnv-example-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.getElementById("cnv-chromosome").value = btn.dataset.chr;
    document.getElementById("cnv-start").value      = btn.dataset.start;
    document.getElementById("cnv-end").value        = btn.dataset.end;
    document.getElementById("cnv-copies").value     = btn.dataset.copies;
    document.getElementById("cnv-type").value       = btn.dataset.type;
    document.querySelectorAll(".cnv-toggle").forEach(b => {
      b.classList.toggle("active", b.dataset.type === btn.dataset.type);
    });
  });
});

// ── Main Submit Button ────────────────────────────────────────────────────────
document.getElementById("cnvBtn").addEventListener("click", async () => {
  const btn  = document.getElementById("cnvBtn");
  const text = document.getElementById("cnvBtnText");
  const spin = document.getElementById("cnvBtnSpinner");

  text.textContent = "Analyzing…";
  spin.classList.remove("hidden");
  btn.disabled = true;

  // Determine active mode
  const activeTab  = document.querySelector(".mode-tab.active");
  const mode       = activeTab ? activeTab.dataset.mode : "gene";
  const cnvType    = document.getElementById("cnv-type").value;
  const copiesVal  = document.getElementById("cnv-copies").value;
  const copyNumber = copiesVal !== "" ? parseInt(copiesVal) : null;

  let url     = "";
  let payload = {};

  if (mode === "rsid") {
    const rsid = document.getElementById("cnv-rsid-input").value.trim().toLowerCase();
    if (!rsid) { showCNVError("Please enter an rsID."); resetCNVBtn(); return; }
    url     = "/api/analyze-cnv-rsid";
    payload = { rsid, cnv_type: cnvType, copy_number: copyNumber };

  } else if (mode === "gene") {
    const gene = document.getElementById("cnv-gene-input").value.trim().toUpperCase();
    if (!gene) { showCNVError("Please enter a gene symbol."); resetCNVBtn(); return; }
    url     = "/api/analyze-cnv-rsid";
    payload = { gene, cnv_type: cnvType, copy_number: copyNumber };

  } else {
    // Coordinate mode
    const chromosome = document.getElementById("cnv-chromosome").value;
    const start      = parseInt(document.getElementById("cnv-start").value);
    const end        = parseInt(document.getElementById("cnv-end").value);
    if (!chromosome || isNaN(start) || isNaN(end)) {
      showCNVError("Please fill in all coordinate fields.");
      resetCNVBtn();
      return;
    }
    if (start >= end) {
      showCNVError("Start position must be less than end position.");
      resetCNVBtn();
      return;
    }
    url     = "/api/analyze-cnv";
    payload = { chromosome, start, end, cnv_type: cnvType, copy_number: copyNumber };
  }

  try {
    const res  = await fetch(url, {
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
      showCNVError(data.error || "Variant not found."); 
    } else { 
      renderCNVResult(data); 
    }
  } catch(err) {
    showCNVError("Network error — is the server running?");
  } finally {
    resetCNVBtn();
  }
});

function resetCNVBtn() {
  document.getElementById("cnvBtnText").textContent = "Analyze CNV";
  document.getElementById("cnvBtnSpinner").classList.add("hidden");
  document.getElementById("cnvBtn").disabled = false;
}

// ── Render CNV Result ─────────────────────────────────────────────────────────
function renderCNVResult(d) {
  document.getElementById("cnvEmptyState").classList.add("hidden");
  const panel = document.getElementById("cnvResultPanel");
  panel.classList.remove("hidden");

  const inp    = d.input;
  const source = inp.source || "coord_lookup";

  // Title
  document.getElementById("cr-title").textContent =
    `Chr${inp.chromosome}: ${inp.start.toLocaleString()}–${inp.end.toLocaleString()} (${inp.cnv_type})`;
  document.getElementById("cr-chr").textContent    = "Chr " + inp.chromosome;
  document.getElementById("cr-region").textContent =
    `${inp.start.toLocaleString()} – ${inp.end.toLocaleString()}`;
  document.getElementById("cr-size").textContent       = d.region_size_bp.toLocaleString() + " bp";
  document.getElementById("cr-size-class").textContent = d.size_class;

  // Source badge
  const srcBadge = document.getElementById("cr-source-badge");
  if (source === "rsid_lookup" && inp.rsid) {
    srcBadge.textContent = inp.rsid;
    srcBadge.classList.remove("hidden");
  } else if (source === "gene_lookup" && inp.gene_symbol) {
    srcBadge.textContent = inp.gene_symbol;
    srcBadge.classList.remove("hidden");
  } else {
    srcBadge.classList.add("hidden");
  }

  // Dosage
  const badge = document.getElementById("cr-dosage-badge");
  badge.textContent = d.dosage_class;
  badge.className   = "impact-badge impact-" + d.dosage_class;
  document.getElementById("cr-dosage-text").textContent = d.dosage_effect;

  // Gene Context & Biological Engine
  if (d.gene_context && d.gene_context.available) {
    document.getElementById("cnvGeneContextContainer").innerHTML = renderGeneContextHTML(d.gene_context, d.research_relevance);
  } else {
    document.getElementById("cnvGeneContextContainer").innerHTML = "";
  }

  if (d.diseases && d.diseases.available) {
    document.getElementById("cnvDiseasesContainer").innerHTML = renderDiseasesHTML(d.diseases);
  } else {
    document.getElementById("cnvDiseasesContainer").innerHTML = "";
  }

  if (d.pathways && d.pathways.available) {
    document.getElementById("cnvPathwaysContainer").innerHTML = renderPathwaysHTML(d.pathways);
  } else {
    document.getElementById("cnvPathwaysContainer").innerHTML = "";
  }

  document.getElementById("cr-clinical-note").textContent = d.clinical_note;

  // Explanations
  document.getElementById("cnvReportSummaryContainer").innerHTML = renderReportSummaryHTML(d);
  document.getElementById("cnvGwasContainer").innerHTML = renderGWASHTML(d);
  document.getElementById("cnvPopulationFrequenciesContainer").innerHTML =
    renderPopulationFrequenciesHTML(d);
  document.getElementById("cnvImpactExplanationContainer").innerHTML =
    renderImpactExplanationHTML(d.impact_explanation, d.dosage_class);
  document.getElementById("cnvSignificanceExplanationContainer").innerHTML =
    renderSignificanceExplanationHTML(d.significance_explanation, d.clinical_significance || "");
  if (d.interpretation_summary) {
    document.getElementById("cnvSummaryContainer").innerHTML = 
      `<div class="gv-card result-card info-card gv-is-3e33972d" >
         <h3 class="gv-card-title">Summary</h3>
         <p  class="gv-is-9e7fd580">${d.interpretation_summary}</p>
       </div>`;
  } else {
    document.getElementById("cnvSummaryContainer").innerHTML = "";
  }

  // References
  document.getElementById("referencesContainer").innerHTML = renderEvidenceHTML(d.evidence);
  document.getElementById("cnvPubmedContainer").innerHTML = renderPubMedHTML(d.pubmed);

  
  document.getElementById("dynamicReportNav")?.classList.remove("hidden");
  document.getElementById("dynamicReportTimestamp")?.classList.remove("hidden");

  panel.scrollIntoView({behavior: "smooth"});
}

function showCNVError(msg) {
  document.getElementById("cnvEmptyState").classList.add("hidden");
  const panel = document.getElementById("cnvResultPanel");
  panel.classList.remove("hidden");
  
  // Clear all previous results
  document.getElementById("cr-title").textContent = "Error";
  document.getElementById("cr-chr").textContent = "—";
  document.getElementById("cr-region").textContent = "—";
  document.getElementById("cr-size").textContent = "—";
  document.getElementById("cr-size-class").textContent = "—";
  document.getElementById("cr-source-badge").classList.add("hidden");
  
  document.getElementById("cr-dosage-badge").textContent = "—";
  document.getElementById("cr-dosage-badge").className = "impact-badge";
  document.getElementById("cr-dosage-text").textContent = "—";
  
  document.getElementById("cr-gene-section").classList.add("hidden");
  document.getElementById("cr-clinical-note").textContent = "—";
  
  document.getElementById("cnvGwasContainer").innerHTML = "";
  document.getElementById("cnvPopulationFrequenciesContainer").innerHTML = "";
  document.getElementById("cnvImpactExplanationContainer").innerHTML = "";
  document.getElementById("cnvSignificanceExplanationContainer").innerHTML = "";
  document.getElementById("cnvSummaryContainer").innerHTML = "";
  document.getElementById("referencesContainer").innerHTML = "";
  document.getElementById("cnvPubmedContainer").innerHTML = "";

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

});
