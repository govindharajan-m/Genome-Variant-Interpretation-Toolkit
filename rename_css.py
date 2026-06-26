import glob
import re

classes_to_rename = [
    'navbar', 'nav-link', 'nav-item', 'container', 'row', 
    'card', 'badge', 'btn', 'table', 'form-control',
    'btn-primary', 'btn-outline', 'btn-mini', 'btn-full',
    'card-title', 'card-subtitle', 'card-arrow'
]

# Update CSS
css_file = 'static/css/style.css'
with open(css_file, 'r', encoding='utf-8') as f:
    css = f.read()

for cls in classes_to_rename:
    css = re.sub(r'\.' + cls + r'\b', '.gv-' + cls, css)

with open(css_file, 'w', encoding='utf-8') as f:
    f.write(css)

# Update HTML templates
for html_file in glob.glob('templates/*.html'):
    with open(html_file, 'r', encoding='utf-8') as f:
        html = f.read()
    
    # Target exact whole words in HTML
    for cls in classes_to_rename:
        html = re.sub(r'(?<=\s)' + cls + r'(?=\s|>|/|"|\')', 'gv-' + cls, html)
        html = re.sub(r'"' + cls + r'(?=\s|>|/|"|\')', '"gv-' + cls, html)
        html = re.sub(r'\'' + cls + r'(?=\s|>|/|"|\')', "'gv-" + cls, html)
    
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html)

print("CSS and HTML files updated successfully.")
