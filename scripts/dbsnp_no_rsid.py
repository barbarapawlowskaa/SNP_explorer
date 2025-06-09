import requests
import time
#Library for reading VCF files
from cyvcf2 import VCF 

#Function to fetch rsID from Ensembl VEP given chromosome position, reference base, and alternate base
def get_rsid_from_ensembl(chrom, pos, ref, alt):
    url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}/{ref}/{alt}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    try:
        r = requests.get(url, headers=headers, timeout=15)
        if not r.ok:
            return None
        result = r.json()
        if result and isinstance(result, list):
            for entry in result:
                #Checking if main identifier starts with 'rs' (meaning it is an rsID)
                if "id" in entry and entry["id"].startswith("rs"):
                    return entry["id"]
                #Checking colocated variants (other variants at the same position)
                colocated = entry.get("colocated_variants", [])
                for coloc in colocated:
                    if "id" in coloc and coloc["id"].startswith("rs"):
                        return coloc["id"]
    except Exception as e:
        print(f"Error: Ensembl request failed for {chrom}:{pos} {alt} → {e}")
    return None

#Function to retrieve clinical data from dbSNP using rsID
def query_dbsnp(rsid):
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
        allele_annotations = primary.get('allele_annotations', [])

        genes = set()
        clinical_significance = "unknown"
        condition = "N/A"
        disease_name = "N/A"

        for annotation in allele_annotations:
            clinical = annotation.get('clinical', [])
            if clinical:
                clinical_significance = clinical[0].get('clinical_significances', ["unknown"])[0]
                condition = clinical[0].get('conditions', [{}])[0].get('preferred_name', "N/A")
                disease_name = clinical[0].get('conditions', [{}])[0].get('name', "N/A")

            #Extracting genes from assembly annotations
            assembly_annotations = annotation.get('assembly_annotation', [])
            for assembly in assembly_annotations:
                for gene in assembly.get('genes', []):
                    if 'gene_symbol' in gene:
                        genes.add(gene['gene_symbol'])
                    if 'locus' in gene:
                        genes.add(gene['locus'])

        return {
            "clinical_significance": clinical_significance,
            "condition": condition,
            "disease_name": disease_name,
            "genes": ";".join(sorted(genes)) if genes else "N/A"
        }
    except Exception as e:
        print(f"Error: dbSNP query failed for {rsid} → {e}")
        return None

#Function to annotate variants from a VCF file and save the result to CSV
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
        #Support for the first ALT allele only
        alt = variant.ALT[0] 

        #rsID from Ensembl
        rsid = get_rsid_from_ensembl(chrom, pos, ref, alt)
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




