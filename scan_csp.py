import os
import re

err_count = 0
tpl_dir = r'd:\genome_variant_platform\genome_variant_platform\Genome-Variant-Interpretation-Toolkit\templates'
for fname in os.listdir(tpl_dir):
    if not fname.endswith('.html'): continue
    path = os.path.join(tpl_dir, fname)
    c = open(path, 'r', encoding='utf-8').read()
    
    scripts = re.findall(r'<script.*?>', c, re.IGNORECASE)
    for s in scripts:
        if 'src=' not in s:
            print(f'{fname}: found inline script block: {s}')
            err_count += 1
            
    events = re.findall(r'(onclick|onchange|onsubmit|onload|oninput)\s*=', c, re.IGNORECASE)
    for e in events:
        print(f'{fname}: found inline event: {e}')
        err_count += 1
        
    if 'document.write' in c:
        print(f'{fname}: found document.write')
        err_count += 1
        
    styles = re.findall(r'style\s*=\s*[\"\']', c, re.IGNORECASE)
    for s in styles:
        print(f'{fname}: found inline style')
        err_count += 1
        
print(f'Total violations found: {err_count}')
