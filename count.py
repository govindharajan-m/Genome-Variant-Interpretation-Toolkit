import re

for fname in ['report.html', 'single_variant.html']:
    path = rf'd:\genome_variant_platform\genome_variant_platform\templates\{fname}'
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        c = f.read()
    styles = re.findall(r'style="([^"]*)"', c)
    print(f'{fname}: {len(styles)}')
