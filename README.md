# SNP_explorer 

SNP Explorer is a command-line tool for analyzing SNP (single nucleotide polymorphism) variants from VCF (Variant Call Format) files. It filters variants based on clinical significance and gene of interest and generates detailed HTML and CSV reports.

## Features

1. Annotating variants using RSIDs or querying Ensembl to resolve missing RSIDs.
2. Querying dbSNP for clinical, frequency and gene annotations.
3. Filtering variants based on presence of RSID, clinical significance, and optional gene input.
4. Generating an HTML summary report and a CSV export of results.

## SNP_explorer installation

```
git clone https://github.com/barbarapawlowskaa/SNP_explorer.git
```
## Dependencies installation

```
cd SNP_explorer
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## SNP_explorer pipeline usage

```
snp_explorer.py [-h] [-c CSV_PATH] [-w] [-wo] [-g GENES] input_vcf report_html
```

where:

```
input_vcf                                input VCF file
report_html                              output HTML report path. The CSV path defaults to the same path with a .csv
                                         extension unless --csv-path (-c) is provided.
-h             --help                    show help message and exit
-c             --csv_path CSV_PATH       optional different output CSV path (default: same as HTML with .csv extension)
-w             --with_rsid               annotate assuming RSIDs are present in VCF
-wo            --without_rsid            annotate assuming no RSIDs in VCF, resolve via Ensembl
-g             --genes GENES             gene filtering
```

## Examples 

```
python3 snp_explorer.py -wo data/sample.vcf reports/report.html -c results/report.csv
```
Annotates assuming there are no RSIDs in the VCF file, stores the HTML report in reports/report.html and saves the CSV in results/report.csv.

```
python3 snp_explorer.py -w -g BRCA1 data/sample.vcf reports/report.html 
```
Annotates assuming RSIDs are present in the VCF file, filters wariants for BRCA1 gene and outputs both the HTML report (as reports/report.html) and CSV report (as reports/report.csv).

```
python3 snp_explorer.py -w -g BRCA1,BRCA2 data/sample.vcf reports/report.html 
```
Annotates assuming RSIDs are present in the VCF file, filters wariants for both BRCA1 and BRCA2 genes and outputs both the HTML report (as reports/report.html) and CSV report (as reports/report.csv).
