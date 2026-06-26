# REMAINING BUGS TRACKER - GenomeVAP V3.0.4

## Status: ZERO KNOWN BUGS
Following the comprehensive V3.0.4 Final Acceptance Testing phase, **there are exactly zero known, reproducible bugs** remaining in the codebase.

### Previously Resolved Vectors (Do Not Reopen)
- **Bootstrap Layout Collapse:** The `gv-row` to `row` and `gv-table` to `table` regressions have been permanently rolled back. Grids now flex properly.
- **Form Data vs JSON Payloads:** Discrepancies between `/api/batch` expecting JSON arrays and the Panel/Cohort engines expecting `multipart/form-data` were handled securely.
- **Python Scope Crash:** `_load_panels` was successfully imported into `variant_engine.py` preventing `NameError` faults.
- **Orphaned Layouts:** All pages missing `<div class="gv-container">` wrappers have been correctly encased to prevent ultra-wide scaling.

### Accepted Anomalies (Not Bugs)
1. **Dark Mode Toggle:** GenomeVAP relies on native CSS `:root` variables overriding standard Bootstrap colors rather than official Bootstrap 5.3+ dark mode switches. This is an intentional design choice to guarantee the phosphor-amber branding cannot be inadvertently washed out by native OS light-mode settings.
2. **Mock Database Bounds:** Certain deep variant queries may gracefully resolve to "Not Found" rather than an error; this represents the functional bounds of the mock internal JSON databases, not an API routing fault.

## System Health
**Ready for deployment.**
