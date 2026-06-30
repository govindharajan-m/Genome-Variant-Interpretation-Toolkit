document.addEventListener("DOMContentLoaded", () => {
let currentPanelData = null;

document.getElementById('panelForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const disease = document.getElementById('diseaseInput').value.trim();
    if (!disease) return;

    document.getElementById('resultsContainer').style.display = 'none';
    document.getElementById('errorContainer').style.display = 'none';
    document.getElementById('loadingIndicator').style.display = 'block';

    try {
        const formData = new FormData();
        formData.append('disease', disease);
        
        const response = await fetch('/api/panel', {
            method: 'POST',
            body: formData
        });
        
        const data = await response.json();
        
        document.getElementById('loadingIndicator').style.display = 'none';
        
        if (!response.ok || data.error) {
            document.getElementById('errorContainer').textContent = data.error || 'An error occurred while generating the panel.';
            document.getElementById('errorContainer').style.display = 'block';
            return;
        }
        
        currentPanelData = data;
        renderPanelResults(data);
        document.getElementById('resultsContainer').style.display = 'block';
        
    } catch (err) {
        document.getElementById('loadingIndicator').style.display = 'none';
        document.getElementById('errorContainer').textContent = 'Network or server error.';
        document.getElementById('errorContainer').style.display = 'block';
    }
});

function toggleVariantView(geneId) {
    const el = document.getElementById('variants-' + geneId);
    if (el.style.display === 'none') {
        el.style.display = 'table-row';
    } else {
        el.style.display = 'none';
    }
}

function renderPanelResults(data) {
    // Summary
    const summary = data.summary;
    document.getElementById('summaryDisease').textContent = summary.disease_name;
    document.getElementById('summaryGeneCount').textContent = summary.genes_evaluated;
    document.getElementById('summaryCritical').textContent = summary.critical_priority_variants;
    document.getElementById('summaryLit').textContent = data.total_literature.toLocaleString() + "+ papers";
    document.getElementById('summaryTopGene').textContent = summary.top_variant;

    // Table
    const tbody = document.getElementById('geneTableBody');
    while(tbody.firstChild) tbody.removeChild(tbody.firstChild);
    
    data.recommended_genes.forEach((g, idx) => {
        const geneId = 'gene-' + idx;
        const tr = document.createElement('tr');
        tr.style.borderBottom = "2px solid var(--border-colour)";
        
        let prioColor = "var(--bg-tertiary)";
        let tierLabel = 'Moderate';
        if (g.panel_score >= 80) { prioColor = "rgba(231, 76, 60, 0.2)"; tierLabel = 'Critical'; }
        else if (g.panel_score >= 50) { prioColor = "rgba(230, 126, 34, 0.2)"; tierLabel = 'High'; }
        else if (g.panel_score >= 25) { prioColor = "rgba(241, 196, 15, 0.2)"; tierLabel = 'Moderate'; }
        
        const td1 = document.createElement('td');
        td1.style.fontWeight = "bold"; td1.style.color = "var(--primary-colour)"; td1.style.fontSize = "1.1rem"; td1.style.verticalAlign = "middle";
        td1.textContent = g.gene + " ";
        const divVar = document.createElement('div');
        divVar.style.cssText = "font-size: 0.85rem; color: #3498db; margin-top: 0.4rem; cursor: pointer; text-decoration: underline; display: inline-block;";
        divVar.onclick = () => toggleVariantView(geneId);
        divVar.innerHTML = `<i class="fas fa-chevron-down"></i> ${g.variant_count} Candidate Variants`;
        td1.appendChild(divVar);
        tr.appendChild(td1);
        
        const td2 = document.createElement('td');
        td2.style.verticalAlign = "middle";
        const spanBadge = document.createElement('span');
        spanBadge.className = "status-badge";
        spanBadge.style.background = prioColor;
        spanBadge.style.border = "1px solid currentColor";
        spanBadge.setAttribute("aria-label", `Priority: ${tierLabel}`);
        spanBadge.textContent = tierLabel;
        td2.appendChild(spanBadge);
        tr.appendChild(td2);
        
        const td3 = document.createElement('td');
        td3.style.fontWeight = "600"; td3.style.fontSize = "1.1rem"; td3.style.verticalAlign = "middle";
        td3.textContent = g.panel_score;
        tr.appendChild(td3);
        
        const td4 = document.createElement('td');
        td4.style.verticalAlign = "middle"; td4.style.color = "var(--text-primary)";
        td4.textContent = g.paper_count.toLocaleString();
        tr.appendChild(td4);
        
        const td5 = document.createElement('td');
        td5.style.verticalAlign = "middle"; td5.style.color = "var(--text-primary)";
        g.reason.split('\n').forEach(r => {
            const spanR = document.createElement('span');
            spanR.style.cssText = "display:block; margin-bottom:0.2rem; font-size:0.9rem;";
            spanR.textContent = r;
            td5.appendChild(spanR);
        });
        tr.appendChild(td5);
        
        tbody.appendChild(tr);
        
        // Nested Variants
        if (g.variants && g.variants.length > 0) {
            const trVariants = document.createElement('tr');
            trVariants.id = 'variants-' + geneId;
            trVariants.style.display = 'none';
            trVariants.style.backgroundColor = 'var(--bg-secondary)';
            
            const tdNested = document.createElement('td');
            tdNested.colSpan = 5;
            tdNested.style.cssText = "padding: 1.5rem; border-bottom: 2px solid var(--border-colour);";
            
            const h5 = document.createElement('h5');
            h5.style.cssText = "margin-bottom: 1rem; color: var(--primary-colour);";
            h5.textContent = `Candidate Variants: ${g.gene}`;
            tdNested.appendChild(h5);
            
            const tableWrap = document.createElement('div');
            tableWrap.className = "table-responsive";
            
            const innerTable = document.createElement('table');
            innerTable.className = "table table-sm";
            innerTable.style.cssText = "background: transparent; margin: 0;";
            innerTable.innerHTML = `
                <thead>
                    <tr  class="gv-is-4cf0ec67">
                        <th>Variant</th>
                        <th>Discovery Score</th>
                        <th>Discovery Tier</th>
                        <th>Discovery Driver</th>
                        <th>Evidence Conf.</th>
                        <th>Research Rel.</th>
                        <th>ACMG Status</th>
                    </tr>
                </thead>
            `;
            const innerTbody = document.createElement('tbody');
            
            g.variants.forEach(v => {
                const trV1 = document.createElement('tr');
                
                const tdV1 = document.createElement('td');
                tdV1.style.fontWeight = "bold"; tdV1.style.color = "var(--primary-colour)";
                tdV1.textContent = v.variant_id;
                trV1.appendChild(tdV1);
                
                const tdV2 = document.createElement('td');
                tdV2.style.fontWeight = "bold";
                tdV2.textContent = v.discovery_score + "/100";
                trV1.appendChild(tdV2);
                
                const tdV3 = document.createElement('td');
                const spanV3 = document.createElement('span');
                spanV3.className = "status-badge badge-strong";
                spanV3.setAttribute("aria-label", `Tier: ${v.discovery_tier}`);
                spanV3.textContent = v.discovery_tier;
                tdV3.appendChild(spanV3);
                trV1.appendChild(tdV3);
                
                const tdV4 = document.createElement('td');
                tdV4.style.fontStyle = "italic";
                tdV4.textContent = v.discovery_driver;
                trV1.appendChild(tdV4);
                
                const tdV5 = document.createElement('td');
                tdV5.textContent = v.evidence_confidence + "/100";
                trV1.appendChild(tdV5);
                
                const tdV6 = document.createElement('td');
                tdV6.textContent = v.research_relevance + "/100";
                trV1.appendChild(tdV6);
                
                const tdV7 = document.createElement('td');
                tdV7.textContent = v.acmg_status;
                trV1.appendChild(tdV7);
                
                innerTbody.appendChild(trV1);
                
                const trV2 = document.createElement('tr');
                const tdV2_1 = document.createElement('td');
                tdV2_1.colSpan = 7;
                tdV2_1.style.cssText = "border-bottom: 1px dashed var(--border-colour); padding-bottom: 1rem; padding-top: 0.5rem;";
                
                const divRea = document.createElement('div');
                divRea.style.cssText = "font-size: 0.85rem; color: var(--text-secondary); margin-bottom: 0.5rem;";
                divRea.innerHTML = "<strong>Why Recommended?</strong><br>";
                v.discovery_reasons.forEach(r => {
                    const tNode = document.createTextNode("✓ " + r);
                    const br = document.createElement('br');
                    divRea.appendChild(tNode);
                    divRea.appendChild(br);
                });
                tdV2_1.appendChild(divRea);
                
                const divFlag = document.createElement('div');
                divFlag.innerHTML = `<strong  class="gv-is-2a56b75d">Research Flags:</strong> `;
                if(v.research_flags && v.research_flags.length > 0) {
                    v.research_flags.forEach(f => {
                        const spanF = document.createElement('span');
                        spanF.className = "badge bg-secondary";
                        spanF.style.marginRight = "0.3rem";
                        spanF.textContent = f;
                        divFlag.appendChild(spanF);
                    });
                } else {
                    const spanNone = document.createElement('span');
                    spanNone.style.cssText = "font-style: italic; color: var(--text-secondary); font-size: 0.8rem;";
                    spanNone.textContent = "None";
                    divFlag.appendChild(spanNone);
                }
                tdV2_1.appendChild(divFlag);
                
                trV2.appendChild(tdV2_1);
                innerTbody.appendChild(trV2);
            });
            
            innerTable.appendChild(innerTbody);
            tableWrap.appendChild(innerTable);
            tdNested.appendChild(tableWrap);
            trVariants.appendChild(tdNested);
            tbody.appendChild(trVariants);
        }
    });
}

// Export CSV
document.getElementById('exportCsvBtn').addEventListener('click', () => {
    if (!currentPanelData) return;
    
    let csvContent = "data:text/csv;charset=utf-8,Gene,Variant,Discovery Score,Discovery Tier,Discovery Driver,Evidence Confidence,Research Relevance,ACMG Status\n";
    currentPanelData.recommended_genes.forEach(g => {
        if (g.variants && g.variants.length > 0) {
            g.variants.forEach(v => {
                let s_gene = sanitizeCSVValue(g.gene);
                let s_vid = sanitizeCSVValue(v.variant_id);
                let s_d_score = sanitizeCSVValue(v.discovery_score);
                let s_d_tier = sanitizeCSVValue(v.discovery_tier);
                let s_d_driver = sanitizeCSVValue(v.discovery_driver);
                let s_ec = sanitizeCSVValue(v.evidence_confidence);
                let s_rr = sanitizeCSVValue(v.research_relevance);
                let s_acmg = sanitizeCSVValue(v.acmg_status);
                csvContent += `${s_gene},${s_vid},${s_d_score},${s_d_tier},"${s_d_driver}",${s_ec},${s_rr},"${s_acmg}"\n`;
            });
        } else {
            let s_gene = sanitizeCSVValue(g.gene);
            csvContent += `${s_gene},None,0,Low,None,0,0,N/A\n`;
        }
    });
const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", `${currentPanelData.disease.replace(/\s+/g, '_')}_Variant_Panel.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
});

// Export JSON
document.getElementById('exportJsonBtn').addEventListener('click', () => {
    if (!currentPanelData) return;
    
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(currentPanelData, null, 2));
    const link = document.createElement("a");
    link.setAttribute("href", dataStr);
    link.setAttribute("download", `${currentPanelData.disease.replace(/\s+/g, '_')}_Variant_Panel.json`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
});

});
