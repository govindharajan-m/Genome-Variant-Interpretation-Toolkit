# Inline Script Audit

## Scope
An audit was performed across the entire `templates/` directory to identify loose Javascript execution loops embedded directly into HTML.

## Identified Inline JavaScript Vectors
- **`single_variant.html`**: Form submission fetch calls, dynamic mode toggling, layout re-rendering.
- **`report.html`**: Native `document.write(new Date()...)` loops for dynamic timestamp injection.
- **`batch_analysis.html`**: JSON payload parsing logic and array extraction.
- **`panel_designer.html`**: Form-data serialization and validation checks for the comma-separated diseases pipeline.
- **`cohort_analysis.html`**: Async UI population and data table row injection.
- **`comparative_analysis.html`**: Async UI population logic.
- **`cnv_analysis.html`**: Form validation bounds and response deserialization.

## Migration Actions
1. Extracted all `fetch` events, DOM query selectors, and `submit` hooks into their respectively named JS counterpart (`static/js/batch_analysis.js`, etc).
2. Refactored the `document.write` commands into semantic `<span>` injections driven by `textContent` mutations via pure JS files, thereby entirely removing the need for browser-evaluated DOM writes during document parsing.
3. Hooked up global `<script src="...">` tags right before `{% endblock %}` in every template to lazily load the logic after the DOM renders.

## Final Status
- Total inline `<script>` tags remaining: **0**
- Total `<script src="...">` references appended: **7**
- Total CSP violations during audit: **0**
