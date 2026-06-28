import os
import re

TEMPLATES_DIR = r"d:\genome_variant_platform\genome_variant_platform\templates"
STATIC_JS_DIR = r"d:\genome_variant_platform\genome_variant_platform\static\js"

if not os.path.exists(STATIC_JS_DIR):
    os.makedirs(STATIC_JS_DIR)

# Exclude base.html since its scripts are already mostly global, or we can handle it manually.
# For others, we extract <script>...</script> into a JS file.
files_to_process = [
    "single_variant.html",
    "report.html"
]

script_regex = re.compile(r"<script>(.*?)</script>", re.DOTALL | re.IGNORECASE)
dynamic_ts_regex = re.compile(r"Last Updated:\s*<script>document\.write\(new Date\(\)\.toISOString\(\)\.replace\('T', ' '\)\.substring\(0, 16\)\);</script>")

for fname in files_to_process:
    path = os.path.join(TEMPLATES_DIR, fname)
    if not os.path.exists(path): continue
    
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    
    # Replace dynamic timestamp script with span
    if dynamic_ts_regex.search(content):
        content = dynamic_ts_regex.sub(r'Last Updated: <span id="dynamicReportTimestampValue"></span>', content)
        ts_code = "\n\n// Timestamp injector\nconst tsEl = document.getElementById('dynamicReportTimestampValue');\nif (tsEl) tsEl.textContent = new Date().toISOString().replace('T', ' ').substring(0, 16);\n"
    else:
        ts_code = ""

    scripts = script_regex.findall(content)
    if not scripts:
        if ts_code:
            js_name = fname.replace(".html", ".js")
            with open(os.path.join(STATIC_JS_DIR, js_name), "w", encoding="utf-8", errors="ignore") as f:
                f.write(ts_code)
            content = content.replace("{% endblock %}", f"<script src=\"{{{{ url_for('static', filename='js/{js_name}') }}}}\"></script>\n{{% endblock %}}")
            with open(path, "w", encoding="utf-8", errors="ignore") as f:
                f.write(content)
        continue

    # Extract all scripts
    combined_js = ""
    for s in scripts:
        combined_js += s.strip() + "\n\n"
    
    combined_js += ts_code

    js_name = fname.replace(".html", ".js")
    
    with open(os.path.join(STATIC_JS_DIR, js_name), "w", encoding="utf-8", errors="ignore") as f:
        f.write(combined_js)
    
    # Replace the FIRST script tag with the src link, remove the others
    # Actually just remove all script tags and put the src before {% endblock %}
    content = script_regex.sub("", content)
    
    # Ensure no empty block scripts that causes issues
    content = content.replace("{% block scripts %}\n\n{% endblock %}", "{% block scripts %}\n<script src=\"{{ url_for('static', filename='js/" + js_name + "') }}\"></script>\n{% endblock %}")
    if "{% block scripts %}" not in content:
        content += "\n{% block scripts %}\n<script src=\"{{ url_for('static', filename='js/" + js_name + "') }}\"></script>\n{% endblock %}\n"
        
    with open(path, "w", encoding="utf-8", errors="ignore") as f:
        f.write(content)
    
    print(f"Processed {fname} -> {js_name}")
