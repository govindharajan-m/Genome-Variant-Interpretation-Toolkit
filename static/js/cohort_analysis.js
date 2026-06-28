let currentCohortData = null;

document.getElementById('cohortForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const rsids = document.getElementById('rsidsInput').value.trim();
    if (!rsids) return;

    document.getElementById('resultsContainer').style.display = 'none';
    document.getElementById('errorContainer').style.display = 'none';
    document.getElementById('loadingIndicator').style.display = 'block';

    try {
        const formData = new FormData();
        formData.append('rsids', rsids);
        
        const response = await fetch('/api/cohort', {
            method: 'POST',
            body: formData
        });
        
        const data = await response.json();
        
        document.getElementById('loadingIndicator').style.display = 'none';
        
        if (!response.ok || data.error) {
            document.getElementById('errorContainer').textContent = data.error || 'An error occurred while analyzing cohort.';
            document.getElementById('errorContainer').style.display = 'block';
            return;
        }
        
        currentCohortData = data;
        renderCohortResults(data);
        document.getElementById('resultsContainer').style.display = 'block';
        
    } catch (err) {
        document.getElementById('loadingIndicator').style.display = 'none';
        document.getElementById('errorContainer').textContent = 'Network or server error.';
        document.getElementById('errorContainer').style.display = 'block';
    }
});

function renderCohortResults(data) {
    const p_analysis = data.pathway_analysis;
    if (p_analysis) {
        document.getElementById('pathwayAnalysisContainer').style.display = 'flex';
        document.getElementById('pathTopGene').textContent = p_analysis.top_gene || 'None';
        document.getElementById('pathTopPathway').textContent = p_analysis.top_pathway || 'None';
        document.getElementById('pathSharedDiseases').textContent = (p_analysis.shared_diseases && p_analysis.shared_diseases.length > 0) ? p_analysis.shared_diseases.join(', ') : 'No shared diseases';
        document.getElementById('pathwayInsightsText').textContent = p_analysis.summary || 'No pathway insights identified.';
    }

    // Summary Updates
    const summary = data.summary;
    document.getElementById('summaryTotal').textContent = summary.total_variants;
    document.getElementById('hlTopGene').textContent = summary.top_gene_representation;
    document.getElementById('summaryCritical').textContent = summary.critical_variants;
    document.getElementById('summaryHigh').textContent = summary.high_variants;
    document.getElementById('summaryAvgEC').textContent = summary.average_evidence_confidence;
    document.getElementById('summaryAvgRR').textContent = summary.average_research_relevance;

    // Table Updates
    const tbody = document.getElementById('cohortTableBody');
    while(tbody.firstChild) tbody.removeChild(tbody.firstChild);
    
    data.variants.forEach((v, idx) => {
        const tr = document.createElement('tr');
        
        let prioColor = "var(--bg-tertiary)";
        if (v.cohort_score >= 80) prioColor = "rgba(231, 76, 60, 0.2)";
        else if (v.cohort_score >= 50) prioColor = "rgba(230, 126, 34, 0.2)";
        else if (v.cohort_score >= 25) prioColor = "rgba(241, 196, 15, 0.2)";
        
        const td1 = document.createElement('td'); td1.style.fontWeight = "bold"; td1.style.verticalAlign = "middle"; td1.textContent = "#" + (idx+1); tr.appendChild(td1);
        const td2 = document.createElement('td'); td2.style.fontWeight = "bold"; td2.style.color = "var(--primary-colour)"; td2.style.fontSize = "1.1rem"; td2.style.verticalAlign = "middle"; td2.textContent = v.variant_id; tr.appendChild(td2);
        const td3 = document.createElement('td'); td3.style.verticalAlign = "middle"; td3.textContent = v.gene; tr.appendChild(td3);
        
        const td4 = document.createElement('td'); td4.style.verticalAlign = "middle";
        const spanBadge = document.createElement('span');
        spanBadge.className = "status-badge";
        spanBadge.style.background = prioColor;
        spanBadge.style.border = "1px solid currentColor";
        spanBadge.setAttribute("aria-label", `Tier: ${v.cohort_tier}`);
        spanBadge.textContent = v.cohort_tier;
        td4.appendChild(spanBadge);
        tr.appendChild(td4);
        
        const td5 = document.createElement('td'); td5.style.fontWeight = "600"; td5.style.fontSize = "1.1rem"; td5.style.verticalAlign = "middle"; td5.textContent = v.cohort_score; tr.appendChild(td5);
        const td6 = document.createElement('td'); td6.style.verticalAlign = "middle"; td6.textContent = v.evidence_confidence; tr.appendChild(td6);
        const td7 = document.createElement('td'); td7.style.verticalAlign = "middle"; td7.textContent = v.research_relevance; tr.appendChild(td7);
        const td8 = document.createElement('td'); td8.style.verticalAlign = "middle"; td8.textContent = v.acmg_status; tr.appendChild(td8);
        const td9 = document.createElement('td'); td9.style.fontStyle = "italic"; td9.style.color = "var(--text-secondary)"; td9.style.verticalAlign = "middle"; td9.textContent = v.discovery_driver; tr.appendChild(td9);
        
        tbody.appendChild(tr);
    });
}

// Export CSV
document.getElementById('exportCsvBtn').addEventListener('click', () => {
    if (!currentCohortData) return;
    
    let csvContent = "data:text/csv;charset=utf-8,Rank,Variant,Gene,Cohort Score,Priority Tier,Evidence Confidence,Research Relevance,ACMG Status
";
    currentCohortData.variants.forEach((v, i) => {
        let s_rank = sanitizeCSVValue(i+1);
        let s_vid = sanitizeCSVValue(v.variant_id);
        let s_gene = sanitizeCSVValue(v.gene);
        let s_score = sanitizeCSVValue(v.cohort_score);
        let s_tier = sanitizeCSVValue(v.cohort_tier);
        let s_ec = sanitizeCSVValue(v.evidence_confidence);
        let s_rr = sanitizeCSVValue(v.research_relevance);
        let s_acmg = sanitizeCSVValue(v.acmg_status);
        csvContent += `${s_rank},${s_vid},${s_gene},${s_score},${s_tier},${s_ec},${s_rr},"${s_acmg}"
`;
    });
const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", `Cohort_Prioritization.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
});

// Export JSON
document.getElementById('exportJsonBtn').addEventListener('click', () => {
    if (!currentCohortData) return;
    
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(currentCohortData, null, 2));
    const link = document.createElement("a");
    link.setAttribute("href", dataStr);
    link.setAttribute("download", `Cohort_Prioritization.json`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
});

