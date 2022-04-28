"""
File toil_processing.py

Contains code for processing RNA-Seq data from Toil project into a working format.
"""

import numpy as np
import pandas as pd


drop_version = lambda x: x[:x.index('.')]
add_primary_isoform = lambda x: x + '-1' if '-' not in x else x
drop_isoform = lambda x: x[:x.index('-')]

trs_annotation = pd.read_csv('../data/raw/ensembl_trs_annot_260422.txt')

trs_annotation = trs_annotation[(~trs_annotation['UniProtKB/Swiss-Prot ID'].isnull()) & (trs_annotation['Transcript type'] == 'protein_coding')]

trs_annotation['uniprot_isoform'] = trs_annotation['UniProtKB isoform ID'].combine_first(trs_annotation['UniProtKB/Swiss-Prot ID'])
trs_annotation['uniprot_isoform'] = trs_annotation['uniprot_isoform'].map(add_primary_isoform)


toil_index = pd.read_csv('../data/raw/toil_expression_trs_260422', sep='\t', usecols=[0,1])

toil_index['sample'] = toil_index['sample'].map(drop_version)
toil_index = toil_index[toil_index['sample'].isin(trs_annotation['Transcript stable ID'])]

trs_uniprot_toil_indices = np.array([-1]+list(toil_index.index)) + 1 #[-1] to account for header row, see next step

toil_data = pd.read_csv('../data/raw/toil_expression_trs_260422', sep='\t', skiprows=lambda x: x not in trs_uniprot_toil_indices)

toil_data = toil_data['sample'].map(drop_version).set_index('sample', drop=True)

toil_data = toil_data.astype('float32') #Switch to float32 from float64 to reduce memory requirements

toil_data = np.around(2 ** toil_data - 0.001) #Convert expression values to TPM
toil_data = toil_data.where(toil_data != -0.0, 0.0)

toil_data['uniprot_isoform'] = toil_data.index.map(trs_annotation.set_index('Transcript stable ID')['uniprot_isoform'])

toil_data = toil_data.groupby('uniprot_id').sum() #Reindex data to UniProt IDs, sum TPM values for transcripts encoding the same protein


toil_data['uniprot_base'] = toil_data.index.map(drop_isoform)

#Leave only genes which have a primary isoform in the dataset
have_primary_isoform = toil_data['uniprot_id'][toil_data.str.contains('-1')].map(drop_isoform).values
toil_data = toil_data[toil_data['uniprot_base'].isin(have_primary_isoform)]

#Leave only genes which have more than one isoform bar primary in the dataset
toil_data = toil_data[toil_data.duplicated('uniprot_base', keep=False)].drop('uniprot_base')

#Convert expression values from TPM to log2(TPM + 1)
toil_data = np.log2(toil_data + 1)


toil_data.to_csv('../data/processed/toil_uniprot_only_splices.csv')