import glob
import re

classes_to_rename = [
    'navbar', 'nav-link', 'nav-item', 'container', 'row', 
    'card', 'badge', 'btn', 'table', 'form-control',
    'btn-primary', 'btn-outline', 'btn-mini', 'btn-full',
    'card-title', 'card-subtitle', 'card-arrow'
]

# Update JS templates
for html_file in glob.glob('static/js/*.js'):
    with open(html_file, 'r', encoding='utf-8') as f:
        html = f.read()
    
    # Target exact whole words in JS
    for cls in classes_to_rename:
        html = re.sub(r'(?<=\s)' + cls + r'(?=\s|>|/|"|\')', 'gv-' + cls, html)
        html = re.sub(r'"' + cls + r'(?=\s|>|/|"|\')', '"gv-' + cls, html)
        html = re.sub(r'\'' + cls + r'(?=\s|>|/|"|\')', "'gv-" + cls, html)
    
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html)

print("JS files updated successfully.")
