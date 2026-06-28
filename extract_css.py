import os
import re
import hashlib

directories = [
    r'd:\genome_variant_platform\genome_variant_platform\templates',
    r'd:\genome_variant_platform\genome_variant_platform\static\js'
]

css_classes = {}
css_file = r'd:\genome_variant_platform\genome_variant_platform\static\css\inline_styles.css'

dynamic_styles = []

def process_file(path):
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        c = f.read()
        
    def tag_replacer(m):
        tag_full = m.group(0)
        # Find style
        style_match = re.search(r'''style=(['"])(.*?)\1''', tag_full)
        if not style_match: return tag_full
        
        style_content = style_match.group(2).strip()
        if not style_content:
            # just remove empty style
            return re.sub(r'''\s*style=(['"]).*?\1''', '', tag_full)
            
        if '{' in style_content or '$' in style_content or '%' in style_content:
            dynamic_styles.append((path, style_content))
            return tag_full
            
        h = hashlib.md5(style_content.encode('utf-8')).hexdigest()[:8]
        cls_name = f"gv-is-{h}"
        css_classes[cls_name] = style_content
        
        # Remove style attribute
        tag_no_style = tag_full[:style_match.start()] + tag_full[style_match.end():]
        tag_no_style = tag_no_style.replace('  ', ' ')
        
        # Add cls_name to class attribute
        class_match = re.search(r'''class=(['"])(.*?)\1''', tag_no_style)
        if class_match:
            existing_classes = class_match.group(2)
            new_classes = existing_classes + ' ' + cls_name
            # Replace old class attribute with new one
            tag_final = tag_no_style[:class_match.start()] + f'class="{new_classes}"' + tag_no_style[class_match.end():]
        else:
            # Inject class just before the closing >
            if tag_no_style.endswith('/>'):
                tag_final = tag_no_style[:-2] + f' class="{cls_name}"/>'
            else:
                tag_final = tag_no_style[:-1] + f' class="{cls_name}">'
                
        return tag_final

    # Regex to match any HTML tag
    # This might match inside jinja or script, but scripts are mostly gone.
    c_new = re.sub(r'<[a-zA-Z0-9_-]+[^>]+>', tag_replacer, c)
    
    if c_new != c:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(c_new)
        print(f"Updated {os.path.basename(path)}")

for d in directories:
    for fname in os.listdir(d):
        if fname.endswith('.html') or fname.endswith('.js'):
            process_file(os.path.join(d, fname))

if css_classes:
    os.makedirs(os.path.dirname(css_file), exist_ok=True)
    with open(css_file, 'w', encoding='utf-8') as f:
        for cls, rule in css_classes.items():
            f.write(f".{cls} {{ {rule} }}\n")
    print(f"Generated {len(css_classes)} static classes in {css_file}")

if dynamic_styles:
    print(f"Found {len(dynamic_styles)} dynamic styles to fix manually:")
    for p, s in dynamic_styles:
        print(f"  {os.path.basename(p)}: {s}")

