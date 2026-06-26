import os
import re

files_to_fix = [
    'templates/panel_designer.html',
    'templates/cohort_analysis.html',
    'templates/comparative_analysis.html'
]

for filepath in files_to_fix:
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        continue
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Add <div class="gv-container"> after {% block content %}
    if '<div class="gv-container">' not in content:
        content = content.replace('{% block content %}', '{% block content %}\n<div class="gv-container">')
        # Close the div before scripts or endblock
        if '<script>' in content:
            content = content.replace('<script>', '</div>\n<script>')
        else:
            content = content.replace('{% endblock %}', '</div>\n{% endblock %}')
            
    # Remove unnecessary inline margins from rows
    content = content.replace('style="margin-top: 2rem; margin-bottom: 2rem;"', '')
    content = content.replace('style="margin-bottom: 2rem;"', '')
    content = content.replace('style="margin-bottom: 2rem; display: none;"', 'style="display: none;"')
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

print("Layouts fixed successfully.")
