document.addEventListener("DOMContentLoaded", () => {
  // ── Tab switching ────────────────────────────────────────────────────────────
  document.querySelectorAll(".batch-tab").forEach(tab => {
    tab.addEventListener("click", () => {
      document.querySelectorAll(".batch-tab").forEach(t => t.classList.remove("active"));
      document.querySelectorAll(".batch-tab-content").forEach(c => c.classList.add("hidden"));
      tab.classList.add("active");
      document.getElementById(tab.dataset.target).classList.remove("hidden");
    });
  });

  // ── File upload handling ─────────────────────────────────────────────────────
  const fileInput = document.getElementById("rsidFile");
  const fileDisplay = document.getElementById("fileNameDisplay");
  fileInput.addEventListener("change", () => {
    if (fileInput.files[0]) {
      fileDisplay.textContent = "📄 " + fileInput.files[0].name;
      fileDisplay.classList.remove("hidden");
    }
  });

  // Drag & drop
  const dropZone = document.getElementById("dropZone");
  dropZone.addEventListener("dragover", e => { e.preventDefault(); dropZone.classList.add("drag-over"); });
  dropZone.addEventListener("dragleave", () => dropZone.classList.remove("drag-over"));
  dropZone.addEventListener("drop", e => {
    e.preventDefault();
    dropZone.classList.remove("drag-over");
    const file = e.dataTransfer.files[0];
    if (file) {
      fileDisplay.textContent = "📄 " + file.name;
      fileDisplay.classList.remove("hidden");
      // Inject into file input via DataTransfer
      const dt = new DataTransfer();
      dt.items.add(file);
      fileInput.files = dt.files;
    }
  });

  // ── Demo loaders ─────────────────────────────────────────────────────────────
  const PATHOGENIC = ["rs334", "rs113993960", "rs1800562", "rs113488022", "rs28897743", "rs80357906"];

  const loadAllBtn = document.getElementById("loadAll");
  loadAllBtn.addEventListener("click", () => {
    document.querySelector(".batch-tab[data-target='textInput']").click();
    const allRsids = loadAllBtn.dataset.rsids;
    document.getElementById("rsidText").value = allRsids;
  });
  document.getElementById("loadPathogenic").addEventListener("click", () => {
    document.querySelector(".batch-tab[data-target='textInput']").click();
    document.getElementById("rsidText").value = PATHOGENIC.join(", ");
  });

  // ── Batch submit ─────────────────────────────────────────────────────────────
  let lastResults = [];

  document.getElementById("batchBtn").addEventListener("click", async () => {
    const btn = document.getElementById("batchBtn");
    const text = document.getElementById("batchBtnText");
    const spin = document.getElementById("batchBtnSpinner");
    text.textContent = "Processing…";
    spin.classList.remove("hidden");
    btn.disabled = true;

    try {
      let response;
      const activeTab = document.querySelector(".batch-tab.active").dataset.target;

      if (activeTab === "fileInput" && fileInput.files[0]) {
        // File upload via FormData
        const formData = new FormData();
        formData.append("rsid_file", fileInput.files[0]);
        response = await fetch("/api/batch", { method: "POST", body: formData });
      } else {
        // Text input
        const text_area = document.getElementById("rsidText").value.trim();
        if (!text_area) { alert("Please enter at least one rsID."); return; }
        const formData = new FormData();
        formData.append("rsid_list", text_area);
        response = await fetch("/api/batch", { method: "POST", body: formData });
      }

      const data = await response.json();
      if (data.error) { alert("⚠ " + data.error); return; }
      lastResults = data.results;
      renderBatchTable(data.results);

    } catch (err) {
      alert("Network error.");
    } finally {
      text.textContent = "Run Batch Analysis";
      spin.classList.add("hidden");
      btn.disabled = false;
    }
  });

  // ── Render table ─────────────────────────────────────────────────────────────
  function renderBatchTable(results) {
    const tbody = document.getElementById("batchTableBody");
    tbody.innerHTML = "";

    // Summary
    const found = results.filter(r => r.found).length;
    const pathogenic = results.filter(r => r.clinical_significance?.toLowerCase().includes("pathogenic")).length;
    const benign = results.filter(r => r.clinical_significance?.toLowerCase() === "benign").length;
    const vus = results.filter(r => r.clinical_significance?.toLowerCase().includes("uncertain")).length;

    document.getElementById("batchSummary").innerHTML = `
    <div class="summary-item"><span class="sum-num">${results.length}</span><span class="sum-lbl">Total</span></div>
    <div class="summary-item found"><span class="sum-num">${found}</span><span class="sum-lbl">Annotated</span></div>
    <div class="summary-item path"><span class="sum-num">${pathogenic}</span><span class="sum-lbl">Pathogenic</span></div>
    <div class="summary-item vus"><span class="sum-num">${vus}</span><span class="sum-lbl">VUS</span></div>
    <div class="summary-item benign"><span class="sum-num">${benign}</span><span class="sum-lbl">Benign</span></div>
  `;

    results.forEach(r => {
      const tr = document.createElement("tr");
      if (!r.found) tr.classList.add("row-not-found");

      const sigClass = getSigClass(r.clinical_significance);
      const impClass = getImpClass(r.impact_level);
      const freq = r.allele_frequency != null ? (r.allele_frequency * 100).toFixed(2) + "%" : "—";

      tr.innerHTML = `
      <td><span class="rsid-tag small">${r.rsid}</span></td>
      <td><strong>${r.gene || "—"}</strong></td>
      <td class="mono">${r.chromosome ? "Chr" + r.chromosome + ":" + r.position?.toLocaleString() : "—"}</td>
      <td class="mono">${r.ref && r.alt ? r.ref + "→" + r.alt : "—"}</td>
      <td><span class="consequence-small">${(r.consequence || "—").replace(/_/g, " ")}</span></td>
      <td><span class="impact-small ${impClass} gv-is-6174c462" title="${r.impact_explanation ? r.impact_explanation.meaning + ' — ' + (r.impact_explanation.interpretation || '') : ''}" >${r.impact_level}</span></td>
      <td>${freq}</td>
      <td><span class="sig-badge ${sigClass} gv-is-6174c462" title="${r.significance_explanation ? r.significance_explanation.definition + ' — ' + (r.significance_explanation.interpretation || '') : ''}" >${r.clinical_significance || "—"}</span></td>
      <td title="${r.condition || ""}">${truncate(r.condition, 30)}</td>
      <td>${r.evidence ? `<span data-color="${r.evidence.evidence_strength === 'High' ? '#27ae60' : r.evidence.evidence_strength === 'Moderate' ? '#f39c12' : '#e74c3c'}" class="fw-semibold dynamic-color-main">${r.evidence.evidence_strength}</span> <span  class="gv-is-1678cde8">(${r.evidence.evidence_score})</span>` : '—'}</td>
      <td>${r.found ? '<a href="/report/' + r.rsid + '" class="gv-btn-mini">Report →</a>' : "—"}</td>
    `;
      tbody.appendChild(tr);
    });

    document.getElementById("batchResults").classList.remove("hidden");
    document.getElementById("batchResults").scrollIntoView({ behavior: "smooth" });
  }

  // ── Utilities ─────────────────────────────────────────────────────────────────
  function getSigClass(sig) {
    if (!sig) return "";
    const s = sig.toLowerCase();
    if (s.includes("pathogenic") && !s.includes("likely")) return "sig-pathogenic";
    if (s.includes("likely_pathogenic")) return "sig-likely-path";
    if (s.includes("uncertain")) return "sig-vus";
    if (s.includes("benign")) return "sig-benign";
    if (s.includes("risk")) return "sig-risk";
    return "";
  }

  function getImpClass(level) {
    const map = {
      HIGH: "impact-HIGH", MODERATE: "impact-MODERATE",
      LOW: "impact-LOW", MODIFIER: "impact-MODIFIER"
    };
    return map[level] || "";
  }

  function truncate(str, n) {
    if (!str) return "—";
    return str.length > n ? str.slice(0, n) + "…" : str;
  }

  // ── Table filter ──────────────────────────────────────────────────────────────
  document.getElementById("tableSearch").addEventListener("input", function () {
    const q = this.value.toLowerCase();
    document.querySelectorAll("#batchTableBody tr").forEach(row => {
      row.style.display = row.textContent.toLowerCase().includes(q) ? "" : "none";
    });
  });

  // ── CSV Download ──────────────────────────────────────────────────────────────
  document.getElementById("downloadCsvBtn").addEventListener("click", async () => {
    if (!lastResults.length) { alert("No results to download."); return; }
    const res = await fetch("/download/csv", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ results: lastResults }),
    });
    const blob = await res.blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = "variant_analysis_report.csv";
    a.click(); URL.revokeObjectURL(url);
  });

  const browseFileBtn = document.getElementById('browseFileBtn');
  if (browseFileBtn) {
      browseFileBtn.addEventListener('click', () => {
          document.getElementById('rsidFile').click();
      });
  }
});
