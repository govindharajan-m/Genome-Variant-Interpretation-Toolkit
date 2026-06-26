# VISUAL BUG TRACKER

## Identified Visual Defects

### 1. Broken Responsive Grid Layouts
- **Page(s):** All pages (e.g., `/`, `/single_variant`, `/panel_designer`, `/comparative_analysis`)
- **Location/Element:** All Bootstrap column containers (`.col-md-6`, `.col-lg-8`, etc.)
- **Defect:** Columns stacked vertically with 100% width instead of aligning side-by-side in a responsive grid. Margins were broken.
- **Root Cause:** A previous automated CSS namespace script blindly renamed the Bootstrap `.row` class to `.gv-row` in the HTML. Because `.gv-row` is not defined as a flex container in `style.css` or Bootstrap, the grid columns lost their flex context and negative margins.
- **Status:** **FIXED**

### 2. Broken Table Styling
- **Page(s):** `/batch_analysis`, `/cohort_analysis`, `/comparative_analysis`
- **Location/Element:** Results tables
- **Defect:** Tables rendered as raw, unstyled HTML tables, losing Bootstrap padding, borders, and structural styling.
- **Root Cause:** The `table` class was renamed to `gv-table`. The custom CSS targets `.data-table` for aesthetics but still relies on the base structural padding of the standard `table` class.
- **Status:** **FIXED**

### 3. Misaligned Badges and Tags
- **Page(s):** Single Variant Report, Batch Analysis
- **Location/Element:** Pathogenicity badges, Impact badges
- **Defect:** Badges lacked proper inline-block padding, font-weight, and border-radius.
- **Root Cause:** The base `.badge` class was renamed to `.gv-badge`. Bootstrap's `.badge` provides the essential structure that custom classes (like `.status-badge`) were implicitly expecting to combine with.
- **Status:** **FIXED**

### 4. Navbar Mobile Toggle & Spacing Collapse
- **Page(s):** Global
- **Location/Element:** Top Navigation Bar (`.navbar-expand-lg`)
- **Defect:** The navbar layout was overly compressed on certain viewports, and the spacing on `.nav-item` elements was lost.
- **Root Cause:** Bootstrap uses `.nav-item` to apply padding/margins to list items inside a navbar. Renaming it to `.gv-nav-item` stripped the spacing.
- **Status:** **FIXED**
