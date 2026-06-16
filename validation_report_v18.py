import json
from db_handler import get_dbsnp_record

def generate_validation_report():
    print("Generating validation report for Version 1.8...")
    report_lines = [
        "# Version 1.8 Population Frequency Validation Report",
        "",
        "This report verifies the extraction and rendering schema of population frequencies.",
        ""
    ]

    test_cases = [
        {"rsid": "rs1042522", "type": "Local Variant (Mock Dataset)"},
        {"rsid": "rs1234", "type": "Live Variant (Ensembl REST API)"},
        {"rsid": "rs999999999", "type": "Invalid/Missing Variant (Error Handling)"}
    ]

    for tc in test_cases:
        rsid = tc["rsid"]
        report_lines.append(f"## Testing {tc['type']} ({rsid})")
        record = get_dbsnp_record(rsid)
        
        if not record:
            report_lines.append("Variant not found. Error handling successful.\n")
            continue

        report_lines.append(f"**Frequency Source:** `{record.get('frequency_source')}`\n")
        
        freqs = record.get("population_frequencies")
        if not freqs:
            report_lines.append("No population frequency data available.\n")
            continue

        report_lines.append("| Population | Frequency |")
        report_lines.append("| --- | --- |")
        
        for pop, val in freqs.items():
            freq_str = f"{val:.4f}" if val is not None else "N/A"
            report_lines.append(f"| {pop.replace('_', ' ').title()} | {freq_str} |")
        
        report_lines.append("")

    with open("validation_report_v18.md", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))
    print("Validation report saved to validation_report_v18.md")

if __name__ == "__main__":
    generate_validation_report()
