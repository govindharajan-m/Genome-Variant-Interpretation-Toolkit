import sys
from pathlib import Path

path = Path(r'd:\genome_variant_platform\genome_variant_platform\static\js\main.js')
c = path.read_text('utf-8')

func = """
function applyDynamicStyles(container) {
    if (!container) return;
    container.querySelectorAll('.dynamic-width-main').forEach(el => {
        el.style.width = el.dataset.width + '%';
    });
    container.querySelectorAll('.dynamic-color-main').forEach(el => {
        el.style.color = el.dataset.color;
    });
    container.querySelectorAll('.dynamic-color-border-main').forEach(el => {
        el.style.color = el.dataset.color;
        el.style.borderColor = el.dataset.color;
    });
}
"""

if "function applyDynamicStyles" not in c:
    path.write_text(func + "\n" + c, 'utf-8')

# Now we must call applyDynamicStyles(document) in single_variant.js, cnv_analysis.js, batch_analysis.js, cohort_analysis.js, comparative_analysis.js, panel_designer.js right after innerHTML is updated!
