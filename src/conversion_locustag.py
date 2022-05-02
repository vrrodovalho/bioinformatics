# -*- coding: utf-8 -*-
"""
@author: vrrodovalho

"""

import os
import sys
import pathlib
import pandas as pd
from os import listdir
from os.path import isfile, join

def construct_mapping_df(input_dir, list_of_file_paths):
    
    df = pd.DataFrame()
    for file in list_of_file_paths:
        sub_df = pd.read_csv(file)
        df = pd.concat([df, sub_df])
    
    return df


def convert(df_2map, df_maping, id_2map=0, counts_2map=1, from_id='', to_id=''):
    df = df_2map.copy()
    mapping = dict(zip(df_maping[from_id], df_maping[to_id]))
    df['id'] = df[0].map(mapping)
    df = df[ ['id', id_2map, counts_2map] ]
    return df

# DIRECTORY SYSTEM
src_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
main_dir = os.path.dirname(src_dir)
root_dir = os.path.dirname(main_dir)
data_dir = pathlib.Path(main_dir) / 'data'
input_dir = pathlib.Path(data_dir) / 'input'
output_dir = pathlib.Path(data_dir) / 'output'
sys.path.insert(0, root_dir)

# Redefine input and output directories
input_dir = input_dir / 'locustag'

# File paths
count_file = input_dir / 'NT.count'
ncbi_tab_file = input_dir / 'proteins_992_249994.csv'


ncbi_files = [f for f in os.listdir(input_dir) if os.path.isfile(
    os.path.join(input_dir, f)) and 'proteins_' in f]
list_of_file_paths = [ os.path.join(input_dir, filename) \
                      for filename in ncbi_files]

df_ncbi = construct_mapping_df(input_dir=input_dir,
                                list_of_file_paths=list_of_file_paths)


# Read files
df_count = pd.read_csv(count_file, sep='\t', header=None)
#df_ncbi = pd.read_csv(ncbi_tab_file)

# convert
df_converted = convert(df_2map=df_count, df_maping=df_ncbi, 
                       from_id='Locus tag', to_id='Protein product')




file = 'J:\\data\\Renan_RNASEQ\\14SM\\gff\\14SM.gff'
with open(file, 'r') as fasta_file:
    lines = fasta_file.readlines()
    
ids = [line.strip().replace('#!genome-build-accession ','').replace('#!genome-build ','').replace('','') \
       for line in lines if line.startswith('#!genome-build')]
ids = list(set(ids))
print(ids)


