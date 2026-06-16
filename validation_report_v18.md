# Version 1.8 Population Frequency Validation Report

This report verifies the extraction and rendering schema of population frequencies.

## Testing Local Variant (Mock Dataset) (rs1042522)
**Frequency Source:** `Local dbSNP Curated Dataset`

| Population | Frequency |
| --- | --- |
| South Asian | 0.7510 |
| European | 0.1012 |
| East Asian | 0.4746 |
| African | 0.3930 |
| American | 0.1731 |
| Global | 0.4000 |

## Testing Live Variant (Ensembl REST API) (rs1234)
**Frequency Source:** `Ensembl REST API (gnomAD / 1000Genomes)`

| Population | Frequency |
| --- | --- |
| Global | 0.1440 |
| South Asian | 0.1967 |
| European | 0.1698 |
| East Asian | 0.1118 |
| African | 0.0851 |
| American | 0.1492 |

## Testing Invalid/Missing Variant (Error Handling) (rs999999999)
**Frequency Source:** `Ensembl REST API (gnomAD / 1000Genomes)`

| Population | Frequency |
| --- | --- |
| Global | 0.0000 |
| South Asian | N/A |
| European | 0.0000 |
| East Asian | N/A |
| African | 0.0000 |
| American | N/A |
