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
