import glob
import re

classes_to_revert = {
    'gv-row': 'row',
    'gv-table': 'table',
    'gv-badge': 'badge',
    'gv-nav-item': 'nav-item'
}

def revert_classes_in_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original_content = content
    for gv_cls, orig_cls in classes_to_revert.items():
        # HTML/JS class replacements
        # We look for the exact gv-class surrounded by spaces, quotes, or tag brackets
        content = re.sub(r'(?<=\s)' + gv_cls + r'(?=\s|>|/|"|\')', orig_cls, content)
        content = re.sub(r'"' + gv_cls + r'(?=\s|>|/|"|\')', '"' + orig_cls, content)
        content = re.sub(r'\'' + gv_cls + r'(?=\s|>|/|"|\')', "'" + orig_cls, content)
        
        # CSS replacements
        content = re.sub(r'\.' + gv_cls + r'\b', '.' + orig_cls, content)

    if content != original_content:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)

# Update HTML templates
for html_file in glob.glob('templates/*.html'):
    revert_classes_in_file(html_file)

# Update JS templates
for js_file in glob.glob('static/js/*.js'):
    revert_classes_in_file(js_file)

# Update CSS
revert_classes_in_file('static/css/style.css')

print("Reverted specific layout classes successfully.")
