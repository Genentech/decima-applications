import pandas as pd
import os, sys

from grelu.data.preprocess import filter_blacklist, filter_chromosomes
from grelu.variant import filter_variants


# See 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/data_tables/dataset_metadata.tsv' for the original metadata
onek1k_cell_type_mapping = {
    'QTD000606': 'B cell',
    'QTD000607': ['memory B cell', 'B cell'],
    'QTD000608': ['naive B cell', 'B cell'],
    'QTD000609': 'classical monocyte',
    'QTD000610': 'non-classical monocyte',
    'QTD000611': 'CD4-positive, alpha-beta T cell',
    'QTD000612': 'CD4-positive, alpha-beta T cell',
    'QTD000613': 'CD4-positive, alpha-beta T cell',
    'QTD000614': 'CD4-positive, alpha-beta T cell',
    'QTD000615': 'CD8-positive, alpha-beta T cell',
    'QTD000616': 'CD8-positive, alpha-beta T cell',
    'QTD000617': 'CD8-positive, alpha-beta T cell',
    'QTD000618': 'hematopoietic stem cell',
    'QTD000619': 'mucosal invariant T cell',
    'QTD000620': 'natural killer cell',
    'QTD000621': 'natural killer cell',
    'QTD000622': 'natural killer cell',
    'QTD000623': 'plasmablast',
    'QTD000624': 'platelet',
    'QTD000625': 'regulatory T cell',
    'QTD000626': 'conventional dendritic cell',
    'QTD000629': 'plasmacytoid dendritic cell',
}

#!wget -r -e robots=off -P /home/karollua/projects/Decima/scborzoi/AKv1/data/eQTL https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000038 
ONEK1K_SUSIE_DIR = '/gstore/data/resbioai/grelu/decima/onek1k/susie/QTS000038'


def load_susie(susie_dir=ONEK1K_SUSIE_DIR):
    susie_df = []
    for ct in os.listdir(susie_dir):
        df = pd.read_table(f'{susie_dir}/{ct}/{ct}.credible_sets.tsv.gz')
        df['chrom'] = [x.split('_')[0] for x in df.variant]
        df['pos'] = [int(x.split('_')[1]) for x in df.variant]
        df['ref'] = [x.split('_')[2] for x in df.variant]
        df['alt'] = [x.split('_')[3] for x in df.variant]
        df['dataset_id'] = ct
        susie_df.append(df)
    
    susie_df = pd.concat(susie_df, axis=0).reset_index(drop=True)
    return susie_df


def filter_susie(susie_df, ad, cell_type_mapping=onek1k_cell_type_mapping):
        
    # Map dataset names
    susie_df = susie_df[susie_df.dataset_id.isin(cell_type_mapping)]
    susie_df['cell_type'] = susie_df.dataset_id.map(cell_type_mapping)

    # Map gene names
    susie_df = susie_df.merge(ad.var.reset_index(names='gene')[['gene_id', 'gene']], on='gene_id', how='left')

    # filter SNVs
    susie_df = filter_variants(susie_df, max_del_len=0, max_insert_len=0, standard_bases=True)

    # filter chromosomes
    susie_df = filter_chromosomes(susie_df, include='autosomesXY')

    # filter blacklist
    susie_df = filter_blacklist(susie_df, genome="hg38", window=100)

    return susie_df[['chrom', 'pos', 'ref', 'alt', 'cell_type', 'gene', 'gene_id', 'rsid', 'region', 'pip', 'pvalue', 'beta']]


# Map cell-types
cell_type_mapping = \
'''model_celltype	eqtl_celltype
B cell	B cell
naive B cell	B cell
memory B cell	memory B cell
CD4-positive, alpha-beta T cell	CD4+ T cell
CD4-positive, alpha-beta T cell	CD4+ CTL cell
CD4-positive, alpha-beta T cell	CD4+ TCM cell
CD4-positive, alpha-beta T cell	CD4+ TEM cell
CD8-positive, alpha-beta T cell	CD8+ T cell
CD8-positive, alpha-beta T cell	CD8+ TCM cell
CD8-positive, alpha-beta T cell	CD8+ TEM cell
regulatory T cell	Treg memory
mucosal invariant T cell	MAIT cell
conventional dendritic cell	dendritic cell
plasmacytoid dendritic cell	plasmacytoid dendritic cell
hematopoietic stem cell	hematopoietic precursor cell
classical monocyte	monocyte
non-classical monocyte	CD16+ monocyte
natural killer cell	NK cell
natural killer cell	CD56+ NK cell
plasmablast	plasmablast
platelet	platelet'''