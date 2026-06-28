import sys
from pathlib import Path

path = Path(r'd:\genome_variant_platform\genome_variant_platform\static\js\main.js')
c = path.read_text('utf-8')

observer_code = """
document.addEventListener('DOMContentLoaded', () => {
    applyDynamicStyles(document.body);
    const observer = new MutationObserver(mutations => {
        mutations.forEach(m => {
            if (m.addedNodes.length) {
                applyDynamicStyles(document.body);
            }
        });
    });
    observer.observe(document.body, { childList: true, subtree: true });
});
"""

if "MutationObserver" not in c:
    path.write_text(c + "\n" + observer_code, 'utf-8')
