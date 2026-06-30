document.addEventListener('DOMContentLoaded', () => {
let currentComparisonData = null;

document.getElementById('compareForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const rsids = document.getElementById('rsidsInput').value.trim();
    if (!rsids) return;

    document.getElementById('resultsContainer').style.display = 'none';
    document.getElementById('errorContainer').style.display = 'none';
    document.getElementById('loadingIndicator').style.display = 'block';

    try {
        const formData = new FormData();
        formData.append('rsids', rsids);
        
        const response = await fetch('/api/compare', {
            method: 'POST',
            body: formData
        });
        
        const data = await response.json();
        
        document.getElementById('loadingIndicator').style.display = 'none';
        
        if (!response.ok || data.error) {
            document.getElementById('errorContainer').textContent = data.error || 'An error occurred while comparing variants.';
            document.getElementById('errorContainer').style.display = 'block';
            return;
        }
        
        currentComparisonData = data;
        renderComparisonResults(data);
        document.getElementById('resultsContainer').style.display = 'block';
        
    } catch (err) {
        document.getElementById('loadingIndicator').style.display = 'none';
        document.getElementById('errorContainer').textContent = 'Network or server error.';
        document.getElementById('errorContainer').style.display = 'block';
    }
});

function createBadge(text, colorClass) {
    const b = document.createElement('span');
    b.className = "badge " + colorClass;
    b.style.display = "block";
    b.style.marginTop = "0.3rem";
    b.textContent = text;
    return b;
}

function renderComparisonResults(data) {
    const summary = data.summary;
    const variants = data.variants;
    const p_analysis = data.pathway_analysis;

    if (p_analysis) {
        document.getElementById('pathwayAnalysisContainer').style.display = 'flex';
        document.getElementById('pathTopGene').textContent = p_analysis.top_gene || 'None';
        document.getElementById('pathTopPathway').textContent = p_analysis.top_pathway || 'None';
        document.getElementById('pathSharedDiseases').textContent = (p_analysis.shared_diseases && p_analysis.shared_diseases.length > 0) ? p_analysis.shared_diseases.join(', ') : 'No shared diseases';
        document.getElementById('pathwayInsightsText').textContent = p_analysis.summary || 'No pathway insights identified.';
    }

    // Update Summary Boxes
    const setSummaryText = (id, arr) => {
        const el = document.getElementById(id);
        el.textContent = (arr && arr.length > 0) ? arr.join(", ") : "None";
    };
    
    setSummaryText('sumTopOverall', summary.top_overall);
    setSummaryText('sumMostClinical', summary.most_clinical);
    setSummaryText('sumMostResearch', summary.most_research);
    setSummaryText('sumMostPgx', summary.most_pgx);
    setSummaryText('sumMostStudied', summary.most_studied);

    // Update Insights
    const ul = document.getElementById('insightsList');
    ul.innerHTML = '';
    if (summary.insights && summary.insights.length > 0) {
        summary.insights.forEach(ins => {
            const li = document.createElement('li');
            li.style.marginBottom = "0.3rem";
            li.textContent = ins;
            ul.appendChild(li);
        });
    } else {
        const li = document.createElement('li');
        li.textContent = "No prominent deterministic insights identified.";
        ul.appendChild(li);
    }

    // Build Table Header
    const thead = document.getElementById('compareThead');
    thead.innerHTML = '';
    const trHead = document.createElement('tr');
    trHead.style.borderBottom = "2px solid var(--border-colour)";
    
    const thMetric = document.createElement('th');
    thMetric.style.width = "20%";
    thMetric.textContent = "Metric";
    trHead.appendChild(thMetric);
    
    variants.forEach(v => {
        const th = document.createElement('th');
        th.style.textAlign = "center";
        th.style.fontSize = "1.2rem";
        th.style.color = "var(--primary-colour)";
        th.textContent = v.variant_id;
        trHead.appendChild(th);
    });
    thead.appendChild(trHead);

    // Build Table Body Rows
    const tbody = document.getElementById('compareTbody');
    tbody.innerHTML = '';
    
    const rowsDef = [
        { label: "Gene", key: "gene" },
        { label: "Clinical Significance", key: "clinical_significance" },
        { label: "Evidence Confidence", key: "evidence_confidence", highlightKey: "Best Clinical Evidence" },
        { label: "Research Relevance", key: "research_relevance", highlightKey: "Best Research Evidence" },
        { label: "Variant Priority", key: "variant_priority", highlightKey: "Highest Overall Priority" },
        { label: "ACMG Status", key: "acmg_status" },
        { label: "Pharmacogenomics", key: "pharmacogenomics", highlightKey: "Strongest Precision Medicine" },
        { label: "Discovery Score", key: "discovery_score" },
        { label: "Cohort Score", key: "cohort_score" },
        { label: "Literature Tier", key: "literature_tier", highlightKey: "Most Studied Variant" },
        { label: "PubMed Count", key: "pubmed_count" },
        { label: "Comparative Strength", key: "_strength" }
    ];

    rowsDef.forEach(r => {
        const tr = document.createElement('tr');
        
        const tdLabel = document.createElement('td');
        tdLabel.style.fontWeight = "bold";
        tdLabel.style.backgroundColor = "var(--bg-secondary)";
        tdLabel.style.color = "var(--text)";
        tdLabel.textContent = r.label;
        tr.appendChild(tdLabel);
        
        variants.forEach(v => {
            const td = document.createElement('td');
            td.style.textAlign = "center";
            td.style.verticalAlign = "middle";
            
            if (r.key === "_strength") {
                const s = v.strength_matrix;
                if (s) {
                    td.innerHTML = `<span  class="gv-is-6623ee52">Strongest: <strong  class="gv-is-b082bd6d">${s.strongest}</strong><br>Weakest: <strong  class="gv-is-cb1de4c7">${s.weakest}</strong></span>`;
                } else {
                    td.textContent = "N/A";
                }
            } else {
                const val = v[r.key];
                td.textContent = val !== null && val !== undefined ? val : "N/A";
                
                // Add badge if this variant won this metric
                if (r.highlightKey && v.badges && v.badges.includes(r.highlightKey)) {
                    td.appendChild(createBadge(r.highlightKey, "bg-warning text-dark"));
                }
            }
            tr.appendChild(td);
        });
        tbody.appendChild(tr);
    });
}

// Export CSV (Vertical Format: Rows=Variants, Cols=Metrics)
document.getElementById('exportCsvBtn').addEventListener('click', () => {
    if (!currentComparisonData) return;
    
    let csvContent = "data:text/csv;charset=utf-8,Variant,Gene,Evidence Confidence,Research Relevance,Variant Priority,Discovery Score,Cohort Score,ACMG Status,PGx Tier\n";
    currentComparisonData.variants.forEach(v => {
        let s_vid = sanitizeCSVValue(v.variant_id);
        let s_gene = sanitizeCSVValue(v.gene);
        let s_ec = sanitizeCSVValue(v.evidence_confidence);
        let s_rr = sanitizeCSVValue(v.research_relevance);
        let s_vp = sanitizeCSVValue(v.variant_priority);
        let s_ds = sanitizeCSVValue(v.discovery_score);
        let s_cs = sanitizeCSVValue(v.cohort_score);
        let s_acmg = sanitizeCSVValue(v.acmg_status);
        let s_pgx = sanitizeCSVValue(v.pharmacogenomics);
        
        csvContent += `${s_vid},${s_gene},${s_ec},${s_rr},${s_vp},${s_ds},${s_cs},"${s_acmg}",${s_pgx}\n`;
    });
    
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", `Comparative_Analysis.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
});

// Export JSON
document.getElementById('exportJsonBtn').addEventListener('click', () => {
    if (!currentComparisonData) return;
    
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(currentComparisonData, null, 2));
    const link = document.createElement("a");
    link.setAttribute("href", dataStr);
    link.setAttribute("download", `Comparative_Analysis.json`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
});

});
