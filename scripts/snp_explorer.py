import argparse
import os
import tempfile
from dbsnp_rsid import annotate_vcf as annotate_with_rsid
from dbsnp_no_rsid import annotate_vcf as annotate_without_rsid
from filter_annotated import filter_variants
from generate_report import generate_report
import pandas as pd

#Main function handling command-line execution of the script
def main():
    #Creating command-line argument parser
    parser = argparse.ArgumentParser(description="SNP Explorer Pipeline")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("report_html", help="Output HTML report path. The CSV path defaults to the same path with a .csv extension unless --csv-path (-c) is provided.")
    parser.add_argument("-c", "--csv_path", help="Optional different output CSV path (default: same as HTML with .csv extension)", default=None)
    #Flags to specify annotation mode: with RSID in VCF, or resolve RSID from Ensembl if missing
    parser.add_argument("-w", "--with_rsid", action="store_true", help="Annotate assuming RSIDs are present in VCF")
    parser.add_argument("-wo", "--without_rsid", action="store_true", help="Annotate assuming no RSIDs in VCF, resolve via Ensembl")
    parser.add_argument("-g", "--genes", help="Gene filtering")
    args = parser.parse_args()

    #Checking file extensions to ensure correct order of positional arguments
    if not args.input_vcf.lower().endswith(".vcf"):
        print("Error: The first positional argument 'input_vcf' should be a .vcf file")
        exit(1)

    if not (args.report_html.lower().endswith(".html") or args.report_html.lower().endswith(".htm")):
        print("Error: The second positional argument 'report_html' should be an .html file")
        exit(1)


    if args.with_rsid == args.without_rsid:
        print("Please specify one of --with_rsid (-w) or --without_rsid (-wo)")
        exit(1)

    #Determining output CSV file path 
    if args.csv_path:
        if os.path.isdir(args.csv_path):
            csv_dir = args.csv_path
            csv_filename = os.path.splitext(os.path.basename(args.report_html))[0] + ".csv"
            final_csv = os.path.join(csv_dir, csv_filename)
        else:
            final_csv = args.csv_path
    else:
        final_csv = os.path.splitext(args.report_html)[0] + ".csv"

    #Creating temporary CSV file to store annotation output
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=True) as tmp_csv:
        if args.with_rsid:
            annotate_with_rsid(args.input_vcf, tmp_csv.name)
        else:
            annotate_without_rsid(args.input_vcf, tmp_csv.name)

        filter_variants(tmp_csv.name, final_csv, gene_filter=args.genes)

    df = pd.read_csv(final_csv)
    sample_name = df.columns[-1] if len(df.columns) > 9 else None
    generate_report(final_csv, args.report_html, sample_name=sample_name)

    print(f"Report generated at: {args.report_html}")
    print(f"CSV report saved at: {final_csv}")

if __name__ == "__main__":
    main()
