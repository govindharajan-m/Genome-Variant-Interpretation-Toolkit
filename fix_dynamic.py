import re
from pathlib import Path

# 1. report.html dynamic styles
report_path = Path(r'd:\genome_variant_platform\genome_variant_platform\templates\report.html')
rc = report_path.read_text('utf-8')

# Replace dynamic widths
def replace_width(m):
    # m.group(0) is the full style="width: {% if ... %}...{% endif %}%;"
    inner = m.group(1)
    return f'data-width="{inner}" class="dynamic-width"'

rc = re.sub(r'style="width:\s*({%.*?%})%;"', replace_width, rc)

# Replace sig_colour
# We have a few patterns:
rc = rc.replace('style="background: {{ report.sig_colour }}22; border: 1px solid {{ report.sig_colour }}; color: {{ report.sig_colour }};"',
                'data-sig-bg="22" data-sig-border="1" data-sig-color="1" class="dynamic-sig-full"')
rc = rc.replace('style="font-size: 1.5rem; font-weight: bold; margin-bottom: 0.5rem; color: {{ report.sig_colour }};"',
                'data-sig-color="1" class="dynamic-sig-color fs-4 fw-bold mb-2"')
rc = rc.replace('style="color: {{ report.sig_colour }}; font-weight: 600;"',
                'data-sig-color="1" class="dynamic-sig-color fw-semibold"')
rc = rc.replace('style="color: {{ report.sig_colour }};"',
                'data-sig-color="1" class="dynamic-sig-color"')
rc = rc.replace('''style="font-weight:600; color:{% if report.research_relevance == 'Very High' %}#8e44ad{% elif report.research_relevance == 'High' %}#2980b9{% elif report.research_relevance == 'Moderate' %}#27ae60{% else %}#7f8c8d{% endif %};"''',
                '''class="fw-semibold color-rel-{{ report.research_relevance | lower | replace(' ', '-') }}"''')
rc = rc.replace('''style="font-weight: bold; color: {% if report.evidence.evidence_strength == 'High' %}#2ecc71{% elif report.evidence.evidence_strength == 'Moderate' %}#f39c12{% else %}#e74c3c{% endif %};"''',
                '''class="fw-bold color-ev-{{ report.evidence.evidence_strength | lower | replace(' ', '-') }}"''')
rc = rc.replace('''style="color: {% if report.impact.impact_level == 'HIGH' %}var(--high){% elif report.impact.impact_level == 'MODERATE' %}var(--moderate){% elif report.impact.impact_level == 'LOW' %}var(--low){% else %}var(--modifier){% endif %};"''',
                '''class="color-impact-{{ report.impact.impact_level | lower }}"''')

# We inject a small script block at the end of report.js or just inline data attributes
# Wait, CSP forbids inline scripts. I will add to report.js
report_path.write_text(rc, 'utf-8')

report_js_path = Path(r'd:\genome_variant_platform\genome_variant_platform\static\js\report.js')
report_js = report_js_path.read_text('utf-8')

report_js += """
document.addEventListener('DOMContentLoaded', () => {
    // We can't access {{ report.sig_colour }} inside JS because Jinja doesn't parse JS files.
    // We should pass it via a meta tag or hidden input.
    const metaSigColour = document.querySelector('meta[name="sig-colour"]');
    const sigColour = metaSigColour ? metaSigColour.getAttribute('content') : '#ffffff';
    
    document.querySelectorAll('.dynamic-width').forEach(el => {
        el.style.width = el.dataset.width + '%';
    });
    document.querySelectorAll('.dynamic-sig-full').forEach(el => {
        el.style.background = sigColour + '22';
        el.style.border = '1px solid ' + sigColour;
        el.style.color = sigColour;
    });
    document.querySelectorAll('.dynamic-sig-color').forEach(el => {
        el.style.color = sigColour;
    });
});
"""
report_js_path.write_text(report_js, 'utf-8')

# We need to inject the meta tag into report.html
rc = rc.replace('<div class="container', '<meta name="sig-colour" content="{{ report.sig_colour }}">\n  <div class="container')
report_path.write_text(rc, 'utf-8')


# 2. main.js and batch_analysis.js
main_path = Path(r'd:\genome_variant_platform\genome_variant_platform\static\js\main.js')
batch_path = Path(r'd:\genome_variant_platform\genome_variant_platform\static\js\batch_analysis.js')

for p in [main_path, batch_path]:
    if not p.exists(): continue
    c = p.read_text('utf-8')
    
    # We replace dynamic styles with classes and we'll append a helper JS function
    # `style="width: ${width}%;"` -> `data-width="${width}" class="dynamic-width-main"`
    c = c.replace('style="width: ${width}%;"', 'data-width="${width}" class="dynamic-width-main"')
    
    # `style="color: ${colour};"` -> `class="main-color-${colour.replace('#', '')}"` -> this doesn't work well if colours are arbitrary.
    # Actually wait, in main.js, what is `colour`?
    # It's returned by getSigColor() which returns "#e74c3c" etc.
    # So we can just use class="main-color-${getSigClass(sig)}" 
    c = c.replace('style="color: ${colour};"', 'data-color="${colour}" class="dynamic-color-main"')
    c = c.replace('style="color: ${color}; border-color: ${color};"', 'data-color="${color}" class="dynamic-color-border-main"')
    c = c.replace("style=\"font-weight:600; color:${relColor};\"", 'data-color="${relColor}" class="fw-semibold dynamic-color-main"')
    c = c.replace("style=\"font-weight: 600; color: ${r.evidence.evidence_strength === 'High' ? '#27ae60' : r.evidence.evidence_strength === 'Moderate' ? '#f39c12' : '#e74c3c'};\"", 'data-color="${r.evidence.evidence_strength === \'High\' ? \'#27ae60\' : r.evidence.evidence_strength === \'Moderate\' ? \'#f39c12\' : \'#e74c3c\'}" class="fw-semibold dynamic-color-main"')
    c = c.replace("style=\"font-weight: 600; color: ${codeColor}; border: 1px solid ${codeColor}; padding: 0.1rem 0.3rem; border-radius: 3px; font-size: 0.8rem;\"", 'data-color="${codeColor}" class="fw-semibold p-1 border rounded fs-6 dynamic-color-border-main"')
    c = re.sub(r'style="padding:\s*0\.75rem\s*0;\s*border-bottom:\s*1px\s*solid\s*var\(--border-2\);\s*\$\{idx\s*===\s*data\.papers\.length\s*-\s*1\s*\?\s*\'border-bottom:\s*none;\'\s*:\s*\'\'\}"', 'class="py-3 border-bottom ${idx === data.papers.length - 1 ? \'border-0\' : \'\'}"', c)
    
    p.write_text(c, 'utf-8')

# We need to run a small observer in main.js to apply data-color to elements injected dynamically.
# Because the HTML strings are injected via innerHTML, we must apply styles right after insertion.
