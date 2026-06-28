

// Timestamp injector
const tsEl = document.getElementById('dynamicReportTimestampValue');
if (tsEl) tsEl.textContent = new Date().toISOString().replace('T', ' ').substring(0, 16);

document.addEventListener('DOMContentLoaded', () => {
    // We can't access {{ report.sig_colour }} inside JS because Jinja doesn't parse JS files.
    // We should pass it via a meta tag or hidden input.
    const metaSigColour = document.querySelector('meta[name="sig-colour"]');
    const sigColour = metaSigColour ? metaSigColour.getAttribute('content') : '#ffffff';
    
    document.querySelectorAll('.dynamic-width').forEach(el => {
        el.style.width = el.dataset.width + '%';
    });
    document.querySelectorAll('.dynamic-sig-full').forEach(el => {
        el.style.background = sigColour + '22';
        el.style.border = '1px solid ' + sigColour;
        el.style.color = sigColour;
    });
    document.querySelectorAll('.dynamic-sig-color').forEach(el => {
        el.style.color = sigColour;
    });

    const printBtn = document.getElementById('printReportBtn');
    if (printBtn) {
        printBtn.addEventListener('click', () => {
            window.print();
        });
    }
});
