import re
from pathlib import Path

# 1. cohort_analysis.html and comparative_analysis.html had `%` which got skipped.
for p in [r'd:\genome_variant_platform\genome_variant_platform\templates\cohort_analysis.html',
          r'd:\genome_variant_platform\genome_variant_platform\templates\comparative_analysis.html']:
    c = Path(p).read_text(encoding='utf-8')
    # just extract the styles as static since they don't have JS/Jinja
    def tag_replacer(m):
        tag_full = m.group(0)
        style_match = re.search(r'''style=(['"])(.*?)\1''', tag_full)
        if not style_match: return tag_full
        
        style_content = style_match.group(2).strip()
        if not style_content: return re.sub(r'''\s*style=(['"]).*?\1''', '', tag_full)
        
        if '{' in style_content or '$' in style_content:
            return tag_full # truly dynamic
        
        import hashlib
        h = hashlib.md5(style_content.encode('utf-8')).hexdigest()[:8]
        cls_name = f"gv-is-{h}"
        
        with open(r'd:\genome_variant_platform\genome_variant_platform\static\css\inline_styles.css', 'a') as f:
            f.write(f".{cls_name} {{ {style_content} }}\n")
            
        tag_no_style = tag_full[:style_match.start()] + tag_full[style_match.end():]
        tag_no_style = tag_no_style.replace('  ', ' ')
        class_match = re.search(r'''class=(['"])(.*?)\1''', tag_no_style)
        if class_match:
            new_classes = class_match.group(2) + ' ' + cls_name
            tag_final = tag_no_style[:class_match.start()] + f'class="{new_classes}"' + tag_no_style[class_match.end():]
        else:
            if tag_no_style.endswith('/>'): tag_final = tag_no_style[:-2] + f' class="{cls_name}"/>'
            else: tag_final = tag_no_style[:-1] + f' class="{cls_name}">'
        return tag_final

    c_new = re.sub(r'<[a-zA-Z0-9_-]+[^>]+>', tag_replacer, c)
    Path(p).write_text(c_new, encoding='utf-8')
    print('Fixed', p)
