"""
File processing_ensembl_annotation.py

Contains code for uniting and processing annotation data for Ensembl transcripts.

Warning: This script is supposed to be run in a notebook from the notebooks directory,
so relative filepaths are specified accordingly.
"""

import datetime

import numpy as np
import pandas as pd

to_bool = lambda x: True if x == 1.0 else False
drop_isoform_na = lambda x: x[:x.index('-')] if x is not np.nan else x
add_primary_isoform = lambda x: x + '-1' if '-' not in x else x


#Use `low_memory=False` to suppress a warning regarding different dtypes in columns
trs_annotation = pd.read_csv('../data/raw/ensembl_annotation_trs_full_20220429.txt', sep='\t', low_memory=False)

trs_annotation = trs_annotation.rename(axis=1, mapper={'Gene stable ID': 'ensembl_gene_id',
                                                       'Transcript stable ID': 'ensembl_trs_id',
                                                       'Protein stable ID': 'ensembl_protein_id',
                                                       'Ensembl Canonical': 'ensembl_is_canonical',
                                                       'Gene name': 'ensembl_gene_name',
                                                       'Transcript type': 'trs_type',
                                                       'Transcript length (including UTRs and CDS)': 'trs_length_bp'})

trs_annotation = trs_annotation.drop_duplicates()

trs_annotation['ensembl_is_canonical'] = trs_annotation['ensembl_is_canonical'].map(to_bool)


trs_uniprot = pd.read_csv('../data/raw/ensembl_annotation_trs_uniprot_20220429.txt', sep='\t')

trs_uniprot = trs_uniprot.rename(axis=1, mapper={'Gene stable ID': 'ensembl_gene_id',
                                                 'Transcript stable ID': 'ensembl_trs_id',
                                                 'UniProtKB isoform ID': 'uniprot_isoform',
                                                 'UniProtKB/Swiss-Prot ID': 'uniprot_base'})

trs_uniprot = trs_uniprot[(~trs_uniprot['uniprot_base'].isnull()) | (~trs_uniprot['uniprot_isoform'].isnull())]

#Drop duplicate rows, then isolate and remove cases where a transcript has multiple UniProt mappings
trs_uniprot = trs_uniprot.drop_duplicates()
trs_uniprot[trs_uniprot['ensembl_trs_id'].duplicated(keep=False)].to_csv('../reports/ensembl_annotation_trs_uniprot_ambiguous.csv', index=False)
trs_uniprot = trs_uniprot.drop_duplicates(subset='ensembl_trs_id', keep=False)

#Ensure there are no cases where a UniProt isoform has no corresponding UniProt base ID
assert trs_uniprot[(~trs_uniprot['uniprot_isoform'].isnull()) & (trs_uniprot['uniprot_base'].isnull())].empty

#Ensure there are no cases where a UniProt isoform has a different base than a corresponding UniProt base ID
trs_uniprot['uniprot_isoform_base'] = trs_uniprot['uniprot_isoform'].map(drop_isoform_na)
assert trs_uniprot[~trs_uniprot['uniprot_isoform_base'].isnull()].query('uniprot_isoform_base != uniprot_base').empty

trs_uniprot['uniprot_isoform'] = trs_uniprot['uniprot_isoform'].combine_first(trs_uniprot['uniprot_base']).map(add_primary_isoform)


trs_annotation['uniprot_base'] = trs_annotation['ensembl_trs_id'].map(trs_uniprot.set_index('ensembl_trs_id')['uniprot_base'])
trs_annotation['uniprot_isoform'] = trs_annotation['ensembl_trs_id'].map(trs_uniprot.set_index('ensembl_trs_id')['uniprot_isoform'])

trs_annotation = trs_annotation[['ensembl_gene_name', 'ensembl_gene_id', 'ensembl_trs_id', 'ensembl_is_canonical',
                                 'trs_type', 'trs_length_bp', 'ensembl_protein_id', 'uniprot_base', 'uniprot_isoform']]

trs_annotation = trs_annotation.sort_values(['ensembl_gene_name', 'uniprot_isoform'], na_position='last')


trs_annotation.to_csv(f'../data/processed/ensembl_annotation_trs_uniprot_{datetime.date.today().strftime("%Y%m%d")}.csv', index=False)