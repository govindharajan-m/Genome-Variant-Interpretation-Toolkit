import os
import re

directories = [
    r'd:\genome_variant_platform\genome_variant_platform\templates',
    r'd:\genome_variant_platform\genome_variant_platform\static\js'
]

unique_styles = set()

for d in directories:
    for fname in os.listdir(d):
        if fname.endswith('.html') or fname.endswith('.js'):
            path = os.path.join(d, fname)
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                c = f.read()
            
            # Match style="..." or style='...'
            styles = re.findall(r'style="([^"]+)"', c)
            styles += re.findall(r"style='([^']+)'", c)
            
            for s in styles:
                unique_styles.add(s)

for idx, s in enumerate(sorted(list(unique_styles))):
    print(f"STYLE_{idx}: {s}")
