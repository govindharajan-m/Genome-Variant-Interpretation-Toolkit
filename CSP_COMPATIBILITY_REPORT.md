# Content Security Policy (CSP) Compatibility Report

## Overview
During the V3.0.4 release, a strict Content Security Policy was enforced which caused critical functionality to fail in the browser because inline Javascript execution and inline CSS styling (`unsafe-inline`) was blocked.

## Root Cause
The `@app.after_request` headers in `app.py` previously defined the following CSP:
`script-src 'self' https://cdn.jsdelivr.net;`
`style-src 'self' 'unsafe-inline' https://cdn.jsdelivr.net https://fonts.googleapis.com;`

Because the `script-src` explicitly omitted `'unsafe-inline'`, all `<script>` tags embedded directly in the Jinja HTML templates were blocked from executing by modern web browsers.

## Resolution Steps
1. **JavaScript Extraction:** We fully decoupled the business logic from the HTML templates. Seven distinct JS payloads were extracted from `single_variant.html`, `report.html`, `panel_designer.html`, `cohort_analysis.html`, `comparative_analysis.html`, `cnv_analysis.html`, and `batch_analysis.html` into independent files inside `/static/js/`.
2. **CSS Abstraction:** We identified and processed over 180 hardcoded inline `style="..."` attributes scattered across the codebase, hashing their properties and assigning them to standardized utility classes generated automatically in `/static/css/inline_styles.css`.
3. **Dynamic Style Fixes:** For styles dynamically mapped to variables (like `% width` values and `sig_colour` HEX codes in `report.html`), we migrated them to `data-` attributes on the DOM elements, then appended a pure CSP-compliant MutationObserver in `main.js` that loops through those attributes and applies the values directly to the CSSOM (`element.style.width = ...`), bypassing CSP constraints entirely.

## Final Result
The CSP in `app.py` has now been strictly enforced:
- `script-src` contains no `unsafe-inline`
- `style-src` contains no `unsafe-inline`

The application is fully functional under strict CSP rules.
