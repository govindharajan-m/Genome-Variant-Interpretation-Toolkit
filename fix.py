import re

path = r'd:\genome_variant_platform\genome_variant_platform\templates\single_variant.html'
with open(path, 'r', encoding='utf-8', errors='ignore') as f:
    content = f.read()

# Fix the dynamic timestamp
content = content.replace(
    "Last Updated: <script>document.write(new Date().toISOString().replace('T', ' ').substring(0, 16));</script>",
    "Last Updated: <span id=\"dynamicReportTimestampValue\"></span>"
)

# Extract and replace the main script
start_idx = content.find('<script>')
end_idx = content.find('</script>', start_idx)

if start_idx != -1 and end_idx != -1:
    content = content[:start_idx] + '<script src=\"{{ url_for(\'static\', filename=\'js/single_variant.js\') }}\"></script>' + content[end_idx+9:]
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)
    print('Replaced main script!')
else:
    print('Main script not found!')
