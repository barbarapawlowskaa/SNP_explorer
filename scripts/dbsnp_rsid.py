from cyvcf2 import VCF
import requests
import time

#Function querying dbSNP API for variant information based on its rsID
def query_dbsnp(rsid):
    #Checking if rsID is valid and starts with "rs"
    if not rsid or not rsid.startswith("rs"):
        return None
    #Removing the 'rs' prefix for use in the API
    snp_id = rsid[2:]
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{snp_id}"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code != 200:
            return None
        data = r.json()
        primary = data.get('primary_snapshot_data', {})
        #Getting allele placement info on the genome
        placements = primary.get('placements_with_allele', [])
        alleles = []
        for placement in placements:
            for allele in placement.get('alleles', []):
                #Get nucleotide sequence of allele
                seq = allele.get('allele', {}).get('spdi', {}).get('inserted_sequence', '')
                if seq and seq not in alleles:
                    alleles.append(seq)

        allele_annotations = primary.get('allele_annotations', [])
        #Dictionary of allele frequencies in populations
        freqs = {}
        clinical_significance = "unknown"
        condition = "N/A"
        disease_name = "N/A"
        genes = set()

        for annotation in allele_annotations:
            #Frequency data
            freq_data = annotation.get('frequency', [])
            for freq_entry in freq_data:
                pop = freq_entry.get('population', {}).get('name', 'unknown')
                af = freq_entry.get('allele_frequency', None)
                if af is not None:
                    freqs[pop] = af

            #Clinical data
            clinical = annotation.get('clinical', [])
            if clinical:
                clinical_significance = clinical[0].get('clinical_significances', ["unknown"])[0]
                condition = clinical[0].get('conditions', [{}])[0].get('preferred_name', "N/A")
                disease_name = clinical[0].get('conditions', [{}])[0].get('name', "N/A")

            #Gene data
            assembly_annotations = annotation.get('assembly_annotation', [])
            for assembly in assembly_annotations:
                genes_data = assembly.get('genes', [])
                for gene in genes_data:
                    if 'locus' in gene:
                        genes.add(gene['locus'])
                    if 'gene_symbol' in gene:
                        genes.add(gene['gene_symbol'])
                
                
        #Returning dictionary with alleles, frequencies, clinical significance, disease info and genes
        return {
            "alleles": alleles,
            "frequencies": freqs,
            "clinical_significance": clinical_significance,
            "condition": condition,
            "disease_name": disease_name,
            "genes": ";".join(sorted(genes)) if genes else "N/A",
        }
    except Exception as e:
        print(f"Errot: dbSNP query failed for {rsid} → {e}")
        return None

#Function to annotate variants from a VCF file and save results to a CSV file
def annotate_vcf(input_vcf, output_file):
    vcf_reader = VCF(input_vcf)
    sample_names = vcf_reader.samples  
    annotated_variants = []

    #Iterate over variants in the VCF file
    for variant in vcf_reader:
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        if not variant.ALT:
            continue
        alt = variant.ALT[0]

        #rsID from Ensembl
        rsid = variant.ID if variant.ID and variant.ID.startswith("rs") else None

        #Clinical info from dbSNP
        dbsnp_info = query_dbsnp(rsid) if rsid else None

        #Collecting sample genotypes
        sample_gts = {}
        genotypes = variant.genotypes
        for sample, gt_info in zip(sample_names, genotypes):
            if gt_info is None or len(gt_info) < 2 or gt_info[0] is None or gt_info[1] is None:
                gt_str = "./."
            else:
                gt_str = f"{gt_info[0]}/{gt_info[1]}"
            sample_gts[sample] = gt_str

        #Saving variant information and annotations
        annotated_variants.append({
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "rsid": rsid if rsid else ".",
            "clinical_significance": dbsnp_info["clinical_significance"] if dbsnp_info else "unknown",
            "condition": dbsnp_info["condition"] if dbsnp_info else "N/A",
            "disease_name": dbsnp_info["disease_name"] if dbsnp_info else "N/A",
            "genes": dbsnp_info["genes"] if dbsnp_info else "N/A",
            "samples": sample_gts
        })

        #Limiting request rate to API (to avoid being blocked)
        time.sleep(0.35)  

    #Saving annotations to output CSV file
    with open(output_file, "w") as f:
        header = [
            "chrom", "pos", "ref", "alt", "rsid",
            "clinical_significance", "condition", "disease_name",
            "genes"
        ] + sample_names
        f.write(",".join(header) + "\n")
        for var in annotated_variants:
            #Row includes variant info plus genotype info for all samples
            row = [str(var[h]) for h in header if h not in sample_names]
            row += [var["samples"][s] for s in sample_names]
            f.write(",".join(row) + "\n")


