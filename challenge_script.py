import vcfpy
import pandas as pd
import requests


VEP_API_URL = 'https://grch37.rest.ensembl.org/variant_effect/vep/human/'

vcf_reader = vcfpy.Reader(open('challenge_data.vcf', 'r'))

annotations = []

for record in vcf_reader:
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt = record.ALT[0] # assuming for now only one alt allele
    # TODO: revisit alts

    # 1. Depth of sequence coverage at site of variation
    depth = record.INFO.get('DP', None)
    allelic_depth = record.INFO.get('AD', None)

    if depth and allelic_depth:
        # 2. Number of reads supporting the variant
        num_reads = allelic_depth[1]

        # 3. Percentage of reads supporting the variant
        #       versus those supporting reference reads
        pct_variant = (num_reads / depth) * 100 if depth > 0 else None
    else:
        num_reads = pct_variant = None
    
    # 4. Query Ensemble VEP API
    vep_params = {
        'variant': f'{chrom}-{pos}-{ref}-{alt}',
        'content-type': 'application/json'
    }

    try:
        response = requests.get(VEP_API_URL, params=vep_params)
    except Exception as e:
        print(e)
        break

    # Parse response
    if response.ok:
        vep_data = response.json()
        gene = vep_data[0]['transcript_consequences'][0]['gene_symbol']
        effect_type = vep_data[0]['transcript_consequences'][0]['consequence_terms'][0]
    else:
        gene = effect_type = None
    
    # Record the annotations
    annotations.append({
        'Chromosome': chrom,
        'Position': pos,
        'Reference': ref,
        'Alternate': alt,
        'Depth': depth,
        'Reads Supporting Variant': num_reads,
        'Percentage Supporting Variant': pct_variant,
        'Gene': gene,
        'Effect Type': effect_type
    })


# Convert annotations to DataFrame
df_annotations = pd.DataFrame(annotations)

# Output to CSV or VCF (depending on preference)
df_annotations.to_csv('annotated_variants.csv', index=False)