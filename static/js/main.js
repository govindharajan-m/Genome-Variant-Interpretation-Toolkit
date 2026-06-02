/**
 * main.js — GenomeVAP / Autoradiograph Terminal
 */

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

// ── References & Evidence HTML Generator ──────────────────────────────────────
function renderReferencesHTML(refs) {
  if (!refs) return "";
  
  let html = `<div class="card result-card" id="referencesCard">
    <h3 class="card-title">Evidence & References</h3>
    <div class="references-container" style="display: flex; flex-direction: column; gap: 1rem;">`;

  const sources = [
    { key: "dbsnp", label: "dbSNP" },
    { key: "clinvar", label: "ClinVar" },
    { key: "pubmed", label: "PubMed" },
    { key: "genecards", label: "GeneCards" }
  ];

  for (const src of sources) {
    const items = refs[src.key] || [];
    html += `<div class="ref-section">
      <div style="font-weight: 600; margin-bottom: 0.25rem; color: var(--text-primary);">${src.label}</div>`;
    
    if (items.length === 0) {
      html += `<div style="color: var(--text-secondary); font-size: 0.875rem;">No reference available</div>`;
    } else {
      html += `<ul style="list-style-type: disc; margin: 0; padding-left: 1.5rem;">`;
      for (const item of items) {
        html += `<li style="margin-bottom: 0.25rem;">
          <a href="${item.url}" target="_blank" class="ext-link" style="color: var(--primary-colour); text-decoration: none;">
            ${item.id}
          </a>
        </li>`;
      }
      html += `</ul>`;
    }
    html += `</div>`;
  }

  html += `</div></div>`;
  return html;
}

