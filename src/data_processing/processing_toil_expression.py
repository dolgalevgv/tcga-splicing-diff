"""
File processing_toil_expression.py

Contains code for processing RNA-Seq data from Toil project into a working format.

Warning: This script is supposed to be run in a notebook from the notebooks directory,
so relative filepaths are specified accordingly.
"""

import datetime

import numpy as np
import pandas as pd


drop_version = lambda x: x[:x.index('.')]
drop_isoform = lambda x: x[:x.index('-')]


trs_annotation = pd.read_csv('../data/processed/ensembl_annotation_trs_uniprot_20220429.csv', low_memory=False)

#Leave only protein-coding transcripts as defined by Ensembl and which also map to a UniProt isoform
trs_annotation = trs_annotation[(trs_annotation['trs_type'] == 'protein_coding') & (~trs_annotation['uniprot_isoform'].isnull())]


toil_index = pd.read_csv('../data/raw/toil_expression_trs_20220426', sep='\t', usecols=[0,1])

toil_index['sample'] = toil_index['sample'].map(drop_version)
toil_index = toil_index[toil_index['sample'].isin(trs_annotation['ensembl_trs_id'])]

toil_trs_uniprot_indices = np.array([-1]+list(toil_index.index)) + 1 #[-1] to account for header row, see next step


toil_data = pd.read_csv('../data/raw/toil_expression_trs_20220426', sep='\t', skiprows=lambda x: x not in toil_trs_uniprot_indices, 
                        index_col=[0], dtype = 'float32', converters = {'sample': str})

toil_data.index = toil_data.index.map(drop_version)

toil_data = np.around(2 ** toil_data - 0.001, decimals=5) #Convert expression values to TPM
toil_data = toil_data.where(toil_data != -0.0, 0.0)

toil_data['uniprot_isoform'] = toil_data.index.map(trs_annotation.set_index('ensembl_trs_id')['uniprot_isoform'])


#Sum TPM values for transcripts encoding the same protein
toil_data = toil_data.groupby('uniprot_isoform').sum()

toil_data['uniprot_base'] = toil_data.index.map(trs_annotation.set_index('uniprot_isoform')['uniprot_base'])

#Leave only genes which have a primary isoform in the dataset
have_primary_isoform = toil_data.index[toil_data.index.str.contains('-1')].map(drop_isoform).unique().values
toil_data = toil_data[toil_data['uniprot_base'].isin(have_primary_isoform)]

#Leave only genes which have more than one isoform bar primary in the dataset
toil_data = toil_data[toil_data.duplicated('uniprot_base', keep=False)].drop('uniprot_base', axis=1)


toil_tcga_data = toil_data[toil_data.columns[toil_data.columns.str.match('TCGA')]]

toil_tcga_data.to_csv(f'../data/processed/toil_tcga_expression_trs_uniprot_splices_only_{datetime.date.today().strftime("%Y%m%d")}.csv')


toil_gtex_data = toil_data[toil_data.columns[toil_data.columns.str.match('GTEX')]]

toil_gtex_data.to_csv(f'../data/processed/toil_gtex_expression_trs_uniprot_splices_only_{datetime.date.today().strftime("%Y%m%d")}.csv')