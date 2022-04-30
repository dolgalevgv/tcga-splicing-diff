"""
File processing_hpa_expression.py

Contains code for processing expression data from Human Protein Atlas.

Warning: This script is supposed to be run in a notebook from the notebooks directory,
so relative filepaths are specified accordingly.
"""

import datetime

import numpy as np
import pandas as pd


trs_annotation = pd.read_csv('../data/processed/ensembl_annotation_trs_uniprot_20220429.csv', low_memory=False)

#Leave only protein-coding transcripts as defined by Ensembl and which also map to a UniProt isoform
trs_annotation = trs_annotation[(trs_annotation['trs_type'] == 'protein_coding') & (~trs_annotation['uniprot_isoform'].isnull())]


hpa_data = pd.read_csv('../data/raw/hpa_cellline_expression_trs_20220408.tsv', sep='\t')

#Map Ensembl trs IDs to UniProt isoforms
hpa_data['uniprot_isoform'] = hpa_data['enstid'].map(trs_annotation.set_index('ensembl_trs_id')['uniprot_isoform'])
hpa_data = hpa_data[~hpa_data['uniprot_isoform'].isnull()]

hpa_data = hpa_data.drop(['ensgid', 'enstid'], axis=1).set_index('uniprot_isoform', drop=True)

hpa_data = hpa_data.groupby('uniprot_isoform').sum()

#Average TPM values for replicates
hpa_data = hpa_data.groupby(hpa_data.columns.str.extract(r'(?<=TPM\.)(.*?)(?=\.)', expand=False), axis=1).mean()

hpa_data = np.log2(hpa_data + 1)


hpa_data.to_csv(f'../data/processed/hpa_cellline_expression_trs_uniprot_only_{datetime.date.today().strftime("%Y%m%d")}.csv')