/**
 * main.js — GenomeVAP Shared JavaScript Utilities
 * ================================================
 * Global helpers loaded on every page via base.html
 */

// ── Auto-uppercase nucleotide inputs ────────────────────────────────────────
document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("input[pattern*='ACGTNacgtn']").forEach(input => {
    input.addEventListener("input", function() {
      const pos = this.selectionStart;
      this.value = this.value.toUpperCase();
      this.setSelectionRange(pos, pos);
    });
  });

  // ── Animate numbers on home page ──────────────────────────────────────────
  document.querySelectorAll(".stat-num").forEach(el => {
    const text = el.textContent;
    // Only animate pure numbers
    const num = parseFloat(text.replace(/[^0-9.]/g, ""));
    if (!isNaN(num) && num > 0 && !text.includes(".")) {
      animateNumber(el, 0, num, 1200, text);
    }
  });
});

/**
 * Animate an element's text content from 0 to a target number.
 * Preserves the surrounding text (e.g. "3.2 Gb", "600M+").
 */
function animateNumber(el, start, end, duration, originalText) {
  // Extract prefix/suffix to preserve "600M+", "20,000+"
  const numStr = end.toLocaleString();
  const prefix = originalText.replace(end.toLocaleString(), "").replace(/[0-9,.]/g,"").trim();
  const suffix = originalText.slice(originalText.lastIndexOf(numStr.slice(-1)) + 1);

  const startTime = performance.now();
  function step(now) {
    const progress = Math.min((now - startTime) / duration, 1);
    const eased    = 1 - Math.pow(1 - progress, 3);
    const current  = Math.floor(start + (end - start) * eased);
    el.textContent = current.toLocaleString() + suffix;
    if (progress < 1) requestAnimationFrame(step);
    else el.textContent = originalText; // ensure final value is exact
  }
  requestAnimationFrame(step);
}

/**
 * Format a genomic position with thousands separator.
 * @param {number} n
 * @returns {string} e.g. "7,676,154"
 */
function formatPos(n) {
  return n?.toLocaleString() ?? "—";
}

/**
 * Derive a CSS class suffix from an impact level string.
 * @param {string} level — "HIGH" | "MODERATE" | "LOW" | "MODIFIER"
 * @returns {string}
 */
function impactClass(level) {
  return "impact-" + (level || "MODIFIER");
}

/**
 * Highlight a search term within a string (returns HTML string).
 * @param {string} text
 * @param {string} query
 * @returns {string}
 */
function highlight(text, query) {
  if (!query || !text) return text || "";
  const re = new RegExp("(" + escapeRe(query) + ")", "gi");
  return text.replace(re, '<mark style="background:rgba(6,182,212,0.25);color:inherit;">$1</mark>');
}

function escapeRe(s) {
  return s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}
