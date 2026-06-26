# FUNCTIONAL BUG TRACKER

## Resolved Bugs

### 1. Panel Designer API Crash (HTTP 500)
- **Endpoint:** `POST /api/panel`
- **Defect:** Server crashed with `NameError: name '_load_panels' is not defined`.
- **Root Cause:** In `variant_engine.py`, the function `_load_panels` was invoked to resolve case-insensitive disease names, but it was not imported from `db_handler.py`.
- **Status:** **FIXED**

### 2. Batch Analysis JSON Payload Crash (HTTP 500)
- **Endpoint:** `POST /api/batch`
- **Defect:** Sending a raw JSON array instead of a JSON object crashed the server.
- **Root Cause:** The route handler assumed the parsed JSON was a dictionary and blindly called `data.get("rsids", [])`, causing an `AttributeError` when a list was provided.
- **Status:** **FIXED** (Payload type checking added).

### 3. Missing Global Layout Constraints
- **Pages:** `/panel_designer`, `/cohort_analysis`, `/comparative_analysis`
- **Defect:** The layout broke out of the centered application boundaries, stretching the full width of the viewport and misaligning with the navigation bar.
- **Root Cause:** The root `<div class="gv-container">` wrapper, which establishes the `max-width` and centering logic, was omitted from these templates.
- **Status:** **FIXED**

### 4. Unnecessary Whitespace
- **Pages:** `/panel_designer`, `/cohort_analysis`, `/comparative_analysis`
- **Defect:** Excessive vertical gaps separated the headers from the interactive forms.
- **Root Cause:** Unnecessary inline styles (`margin-top: 2rem; margin-bottom: 2rem;`) were arbitrarily injected onto the container rows.
- **Status:** **FIXED**

## Action Plan Status
All functional logic and layout constraints have been audited and repaired. The platform is ready for production scaling.
