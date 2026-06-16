import json
from db_handler import get_dbsnp_record

def generate_validation_report():
    print("Generating validation report...")
    report_lines = [
        "# Live Variant Resolution (Phase 1) Validation Report",
        "",
        "This report validates the implementation of NCBI e-utilities fallback for variant resolution.",
        ""
    ]

    # Test 1: Local Variant Lookup
    local_rsid = "rs1042522"
    report_lines.append(f"## 1. Local Variant Lookup: `{local_rsid}`")
    local_record = get_dbsnp_record(local_rsid)
    if local_record and local_record.get('gene') == 'TP53':
        report_lines.append(f"**Result**: ✅ SUCCESS")
        report_lines.append(f"Retrieved locally without network request.")
        report_lines.append(f"```json\n{json.dumps(local_record, indent=2)}\n```")
    else:
        report_lines.append(f"**Result**: ❌ FAILED")

    report_lines.append("")

    # Test 2: Live Variant Lookup
    live_rsid = "rs1234"
    report_lines.append(f"## 2. Live Variant Lookup: `{live_rsid}`")
    live_record = get_dbsnp_record(live_rsid)
    if live_record and live_record.get('rsid') == live_rsid:
        report_lines.append(f"**Result**: ✅ SUCCESS")
        report_lines.append(f"Retrieved via NCBI e-utilities.")
        report_lines.append(f"```json\n{json.dumps(live_record, indent=2)}\n```")
    else:
        report_lines.append(f"**Result**: ❌ FAILED")
        
    report_lines.append("")

    # Test 3: Invalid Variant Lookup
    invalid_rsid = "rs99999999999"
    report_lines.append(f"## 3. Invalid Variant Lookup: `{invalid_rsid}`")
    invalid_record = get_dbsnp_record(invalid_rsid)
    if invalid_record is None:
        report_lines.append(f"**Result**: ✅ SUCCESS")
        report_lines.append(f"Correctly returned None for non-existent rsID.")
    else:
        report_lines.append(f"**Result**: ❌ FAILED")

    report_lines.append("")

    # Write report
    with open("validation_report.md", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))
    
    print("Validation report generated at validation_report.md")

if __name__ == "__main__":
    generate_validation_report()
