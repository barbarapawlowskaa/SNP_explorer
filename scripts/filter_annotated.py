
#Function filtering CSV file of variants based on clinical significance and gene names from RSID column
def filter_variants(input_csv, output_csv, gene_filter=None):
    import pandas as pd

    df = pd.read_csv(input_csv)
    filtered = df[(df["rsid"] != ".") & (df["clinical_significance"] != "unknown")]

    #If gene_filter is provided, further filter variants to those whose "genes" column contains any of the gene names specified in gene_filter
    if gene_filter:
        genes = [g.strip().lower() for g in gene_filter.split(",")]
        filtered = filtered[
            filtered["genes"].str.lower().apply(lambda x: any(g in x for g in genes))
        ]

    #Choosing the columns for CSV file
    base_cols = ["chrom", "pos", "rsid", "clinical_significance", "genes"]
    sample_cols = [col for col in filtered.columns if col not in base_cols and col not in ["ref", "alt", "condition", "disease_name"]]
    if sample_cols:
        base_cols.append(sample_cols[0])


    filtered = filtered[base_cols]
    filtered.to_csv(output_csv, index=False)
