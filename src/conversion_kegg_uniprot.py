# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 18:13:55 2021

@author: Cliente

This script is to convert KEGG and Uniprot ids.

"""

import urllib
import os
import sys
import re
import pathlib

def retrieve_uniprot_2_kegg(format1='hsa', format2='uniprot', mode=1):
    '''
    

    Parameters
    ----------
    format1 : TYPE, optional
        DESCRIPTION. The default is 'hsa'.
    format2 : TYPE, optional
        DESCRIPTION. The default is 'uniprot'.
    mode : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    print("From KEGG, retrieving map for FORMAT CONVERSION...")
    url = 'http://rest.kegg.jp/conv/{}/{}'.format(format1, format2)
    with urllib.request.urlopen(url) as f:
        format_map = f.read().decode('utf-8').splitlines()
    
    print("Parsing format conversion information...")
    pattern = '[\S]+:([\S]+)\t{}:([\d]+)'.format(format1)
    conv1 = {}
    conv2 = {}
    for string in format_map:
        match_obj = re.match(pattern, string)
        (key, value) = (match_obj.group(1), match_obj.group(2)) 
        if key in conv1:
            conv1[key].append(value)
        else:
            conv1[key] = [value]
        if value in conv2:
            conv2[value].append(key)
        else:
            conv2[value] = [key]
    if mode == 1:
        return conv1
    elif mode == 2:
        return conv2
    else:
        print("Choose mode 1 or 2.")
        return None

##############################################################################

# DIRECTORY SYSTEM
src_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
main_dir = os.path.dirname(src_dir)
root_dir = os.path.dirname(main_dir)
data_dir = pathlib.Path(main_dir) / 'data'
input_dir = pathlib.Path(data_dir) / 'input'
output_dir = pathlib.Path(data_dir) / 'output'
sys.path.insert(0, root_dir)

# Redefine input and output directories
input_dir = input_dir / 'fasta_filters'
# output_dir = output_dir / 'fasta_filters'

# Get human uniprot IDs that are in kegg
kegg2uniprot = retrieve_uniprot_2_kegg(format1='hsa', format2='uniprot', 
                                       mode=2)
all_uniprot_ids = []
for kegg_id in kegg2uniprot:
    uniprot_ids = kegg2uniprot[kegg_id]
    all_uniprot_ids.extend(uniprot_ids)
    
all_uniprot_ids = list(set(all_uniprot_ids))

output_file = output_dir / 'human_proteines_list.txt'
with open(output_file, 'w') as f:
    for item in all_uniprot_ids:
        f.write("%s\n" % item)
        
        