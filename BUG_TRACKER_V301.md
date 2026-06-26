# BUG TRACKER V301

## Resolved Bugs
1. **[Critical] `/api/compare` Crash:** `generate_variant_comparison` was not imported in `app.py`. Fixed.
2. **[Critical] Variant Engine Population Crash:** Renamed `eur` references to `cau` inside the variant engine to fix `NameError`.
3. **[High] Bootstrap CSS Collisions:** Bootstrap loaded after the custom `style.css` in `base.html`, overriding layout. Re-ordered correctly and implemented a `gv-` namespace for overlapping classes (`navbar`, `btn`, `card`, etc.).
4. **[High] Render Deployment Failure:** Flask `host` was bound to `127.0.0.1`. Fixed to `0.0.0.0`.
5. **[Medium] Invalid CSS Syntax:** Unresolved pseudo-element (`::bottom`) removed.
6. **[Low] UI Rendering Artifacts:** Fixed `Caucasian/Caucasian` string duplicate typo.
7. **[Low] Undefined CSS Variables:** Fixed references to `--border-1`, `--text-1`, `--bg-elevated`, and `--accent-cyan`.

## Remaining Bugs
- None verified. The platform is stable and all known structural and runtime regressions have been mitigated.

## Risk Assessment
- **Low Risk:** The namespacing of the CSS classes comprehensively insulates GenomeVAP from third-party CSS library updates. 
- **Low Risk:** The codebase relies heavily on local JSON mock databases (`db_handler.py`). True network failures are not fully observable in the current state unless integrated with live APIs.

## Technical Debt
- **Frontend Modularity:** JavaScript logic is heavily consolidated into a single `main.js` file (~800 lines). Future development should isolate component logic (e.g., separating `cohort` logic from `single_variant` logic).
- **CSS Modularity:** `style.css` is massive (~1000 lines). Adopting SCSS/SASS and splitting it into components would improve maintainability.

## Future Recommendations
- Implement comprehensive automated unit testing using `pytest` to prevent `NameError` regressions for undeclared variables or missing imports.
- Re-integrate live API calls natively to transition out of the mock data environment when required.
