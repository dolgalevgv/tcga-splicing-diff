"""
File processing_depmap_expression.py

Contains code for processing expression data from DepMap.

Warning: This script is supposed to be run in a notebook from the notebooks directory,
so relative filepaths are specified accordingly.
"""

import datetime

import numpy as np
import pandas as pd


extract_gene_name = lambda x: x[:x.index(' ')]
extract_trs = lambda x: x[x.index(' ')+2:-1]


trs_annotation = pd.read_csv('../data/processed/ensembl_annotation_trs_uniprot_20220429.csv', low_memory=False)

#Leave only protein-coding transcripts as defined by Ensembl and which also map to a UniProt isoform
trs_annotation = trs_annotation[(trs_annotation['trs_type'] == 'protein_coding') & (~trs_annotation['uniprot_isoform'].isnull())]


depmap_info = pd.read_csv('../data/raw/depmap_22q1_cellline_info_20220413.csv')

depmap_info = depmap_info.set_index('DepMap_ID')['stripped_cell_line_name']


depmap_data = pd.read_csv('../data/raw/depmap_22q1_expression_trs_20220413.csv')

depmap_data = depmap_data.set_index('Unnamed: 0', drop=True).transpose()

depmap_data.columns.name = None
#Replace cell lines' DepMap ID with full names and sort for clearance
depmap_data.columns = depmap_data.columns.map(depmap_info)
depmap_data = depmap_data[sorted(depmap_data.columns)]

#Drop non-ENST entries and separate index into gene name and Ensembl transcript ID
depmap_data = depmap_data[depmap_data.index.str.contains('\(ENST')]
depmap_data['gene_name'] = depmap_data.index.map(extract_gene_name)
depmap_data['ensembl_trs_id'] = depmap_data.index.map(extract_trs)

#Map Ensembl trs IDs to UniProt isoforms
depmap_data['uniprot_isoform'] = depmap_data['ensembl_trs_id'].map(trs_annotation.set_index('ensembl_trs_id')['uniprot_isoform'])
depmap_data = depmap_data[~depmap_data['uniprot_isoform'].isnull()]

depmap_data = depmap_data.drop(['gene_name', 'ensembl_trs_id'], axis=1).set_index('uniprot_isoform', drop=True)

#Convert expression values from log2(TPM +1) to TPM
depmap_data = 2 ** depmap_data - 1

#Sum TPM values for transcripts encoding the same protein
depmap_data = depmap_data.groupby('uniprot_isoform').sum()

#Convert expresion values back to log2(TPM + 1) 
depmap_data = np.log2(depmap_data + 1)


depmap_data.to_csv(f'../data/processed/depmap_22q1_expression_trs_uniprot_only_{datetime.date.today().strftime("%Y%m%d")}.csv')