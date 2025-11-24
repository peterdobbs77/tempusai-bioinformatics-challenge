import vcfpy
import pandas as pd
import requests
import json


VEP_API_URL = 'https://grch37.rest.ensembl.org/vep/human/'

vcf_reader = vcfpy.Reader(open('challenge_data.vcf', 'r'))

annotations = []

for record in vcf_reader:
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt = record.ALT[0] # assuming for now only one alt allele
    # TODO: revisit alts here and below for `alt_depth`
    # alt = ",".join(str(_alt) for _alt in record.ALT)

    # 1. "DP": Depth of sequence coverage at site of variation
    depth = record.INFO.get('DP', None)

    # 3. Percentage of reads supporting the variant
    #       versus those supporting reference reads
    ref_depth = record.INFO.get('RO', None)
    alt_depth = record.INFO.get('AO', None)
    if ref_depth and alt_depth:
        total_depth = ref_depth + alt_depth[0] # NOTE: This number should be same as `depth` above
        pct_variant = (alt_depth[0] / total_depth) * 100 if total_depth > 0 else 0
    else:
        pct_variant = None
    
    # 2. Number of reads supporting the variant
    if depth and pct_variant:
        variant_count = depth*pct_variant
    else:
        variant_coun = None
    
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
        # with open("script_response.txt", "a") as fp:
        #     fp.write("\n")
        #     fp.write(json.dump(response.json()))
        vep_data = response.json()

        if 'transcript_consequences' in vep_data[0]:
            gene = vep_data[0]['transcript_consequences'][0]['gene_symbol']
            effect_type = vep_data[0]['transcript_consequences'][0]['consequence_terms'][0]
        else:
            gene = effect_type = None

        most_severe_consequence = vep_data[0]['most_severe_consequence']
    else:
        gene = effect_type = most_severe_consequence = None
    
    # 5. Minor Allele Frequency (MAF)
    allele_frequency = record.INFO.get('AF', None)

    # Record the annotations
    annotations.append({
        'Chromosome': chrom,
        'Position': pos,
        'Reference': ref,
        'Alternate': alt.value,
        'Alternate Type': alt.type,
        'Depth': depth,
        'Allele Frequency': allele_frequency,
        'Ref Depth': ref_depth,
        'Alt Depth': alt_depth[0],
        'Percentage Supporting Variant': pct_variant,
        'Count Supporting Variant': variant_count,
        'Gene': gene,
        'Effect Type': effect_type,
        'Most Severe Consequence': most_severe_consequence
    })


# Convert annotations to DataFrame
df_annotations = pd.DataFrame(annotations)

# Output to CSV or VCF (depending on preference)
df_annotations.to_csv('annotated_variants.csv', index=False)