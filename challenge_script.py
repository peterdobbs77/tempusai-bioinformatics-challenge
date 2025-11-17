import vcfpy
import pandas as pd
import requests


VEP_API_URL = 'https://grch37.rest.ensembl.org/vep/human/'

vcf_reader = vcfpy.Reader(open('challenge_data.vcf', 'r'))

annotations = []

for record in vcf_reader:
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt = record.ALT[0] # assuming for now only one alt allele
    # TODO: revisit alts
    # alt = ",".join(str(_alt) for _alt in record.ALT)

    # 1. "DP": Depth of sequence coverage at site of variation
    depth = record.INFO.get('DP', None)
    # 2. "AD": Number of reads supporting the variant
    allelic_depth = record.INFO.get('AD', None)

    if depth and allelic_depth:
        # 3. Percentage of reads supporting the variant
        #       versus those supporting reference reads
        ref_depth = allelic_depth[0]
        alt_depth = allelic_depth[1] if len(allelic_depth) > 0 else 0
        total_depth = ref_depth + alt_depth

        pct_variant = (alt_depth / total_depth) * 100 if total_depth > 0 else 0
    else:
        ref_depth = alt_depth = pct_variant = None
    
    # 4. Query Ensemble VEP API
    vep_extension = f"region/{chrom}:{pos}/{alt.value}"
    vep_params = {
        'content-type': 'application/json'
    }
    try:
        response = requests.get(VEP_API_URL+vep_extension, params=vep_params)
    except Exception as e:
        print(e)
        break

    # Parse response
    if response.ok:
        with open("script_response.txt", "a") as fp:
            fp.write("\n")
            fp.write(response)
        vep_data = response.json()
        gene = vep_data[0]['transcript_consequences'][0]['gene_symbol']
        effect_type = vep_data[0]['transcript_consequences'][0]['consequence_terms'][0]
        most_severe_consequence = vep_data[0]['most_severe_consequence']

        # gene = vep_data.get('genes', [{'id': 'N/A'}])[0]['id'] if vep_data else 'N/A'
        # variant_type = vep_data.get('variant_class', 'N/A') if vep_data else 'N/A'
        # effect = vep_data.get('most_severe_consequence', 'N/A') if vep_data else 'N/A'

        # 5. Minor Allele Frequency (MAF)
        # TODO
    else:
        gene = effect_type = most_severe_consequence = None
    
    # Record the annotations
    annotations.append({
        'Chromosome': chrom,
        'Position': pos,
        'Reference': ref,
        'Alternate': alt.value,
        'Alternate Type': alt.type,
        'Depth': depth,
        'Ref Depth': ref_depth,
        'Alt Depth': alt_depth,
        'Percentage Supporting Variant': pct_variant,
        'Gene': gene,
        'Effect Type': effect_type,
        'Most Severe Consequence': most_severe_consequence
    })


# Convert annotations to DataFrame
df_annotations = pd.DataFrame(annotations)

# Output to CSV or VCF (depending on preference)
df_annotations.to_csv('annotated_variants.csv', index=False)