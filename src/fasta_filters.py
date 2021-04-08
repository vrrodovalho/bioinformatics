# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:48:34 2019

@author: vrrodovalho

This script XXXXXXXXXXXXX

"""

import os
import sys
import pathlib


def fasta_parser(fasta_file, forbidden_lines=['',' '], verbose=True):
    '''
    Parses a fasta file into a dictionary, accounting for duplicated ids.

    Parameters
    ----------
    fasta_file : path.
        path of the fasta file.
    forbidden_lines : STR, optional
        Lines that should be ignored. The default is ['',' '].
    verbose : BOOL, optional
        Describe number of sequences if true. The default is True.

    Returns
    -------
    final : Dictionary
        A dictionary in which fasta_ids are the keys and the sequences
        are the values.

    '''

    sequences = {}
    duplicated_ids = {}
        
    with open(fasta_file, 'r') as input_file:
        for line in input_file:
            if line.startswith(">"):
                fasta_id = line[1:].strip()
                if fasta_id in sequences:
                    duplicated_ids[fasta_id] = ""
                    duplicated = True
                else:
                    sequences[fasta_id] = ""
                    duplicated = False
            elif line in forbidden_lines:
                pass
            elif duplicated:
                duplicated_ids[fasta_id] += line.strip()
            else:
                sequences[fasta_id] += line.strip()
    if duplicated_ids:
        final = (sequences, duplicated_ids)
        if verbose:
            print("\nGot {} unique fasta ids!".format(len(sequences))) 
            print("Got {} duplicated fasta ids!".format(len(duplicated_ids)))
    else:
        final = sequences
        if verbose:
            print("\nGot {} unique fasta ids!".format(len(sequences))) 
    return final

def list_parser(list_file, forbidden_lines=['',' '], verbose=True):
    '''
    Parse a file containing a list of ids.

    Parameters
    ----------
    list_file : PATH
        The path of a list file.
    forbidden_lines : LIST, optional
        Lines that should be ignored. The default is ['',' '].
    verbose : BOOL, optional
        Whether to print or not informative strings about the filtering 
        process. The default is True.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    ids = []
    duplicated_ids = []
        
    with open(list_file, 'r') as input_file:
        for line in input_file:
            line = line.strip()
            if line not in forbidden_lines:
                if line not in ids:
                    ids.append(line)
                else:
                    duplicated_ids.append(line)
    len_ids = len(ids)
    len_dup = len(duplicated_ids)
    if duplicated_ids:
        if verbose:
            print("\nGot {} unique ids and {} duplicated ids".format(
                   len_ids, len_dup))
        return (ids, duplicated_ids)
    else:
        if verbose:
            print("\nGot {} unique ids in the list.".format(len_ids))
        return ids

def filter_sequences(sequences, 
                     min_seq_len=35,
                     forbidden=['B','J','O','U','X','Z'],
                     id_filters=[],
                     filter_by={'seq_size':False,
                                'seq_char':False,
                                'id_list_in':False,
                                'id_list_out':False},
                     explain=True):
    '''
    Filter a dictionary containing fasta ids and sequences, based on
    sequence size, sequence alphabet, and allowed of forbidden ids.

    Parameters
    ----------
    sequences : DICT
        Fasta dictionary.
    min_seq_len : INT, optional
        Minimum size accepted for the sequences. The default is 35.
    forbidden : LIST, optional
        List of forbidden characters for the sequences. 
        The default is ['B','J','O','U','X','Z'].
    id_filters : LIST, optional
        List of strings to be applied in fasta_id filters. 
        The default is [].
    filter_by : DICT, optional
        Dictionary specifying which filters should be applied. 
        The default is {'seq_size':False,
                        'seq_char':False,
                        'id_list_in':False,
                        'id_list_out':False}.
        seq_size    Refers to filter by size. 
                    If True, min_seq_len should be specified.
        seq_char    Refers to filter out sequences containing forbidden 
                    characters. If True, forbidden should be specified.
        id_list_in  Filter keeping only fasta_ids containing or matching
                    strings in id_list_in 
        id_list_out Filter taking out only fasta_ids containing or matching
                    strings in id_list_out
        Attention: id_list_out is prioritary over id_list_in.
    explain : BOOL, optional
        Whether to print or not informative strings about the filtering 
        process. The default is True.

    Returns
    -------
    DICT
        A dictionary if no filters were applied.
        Otherwise, two dictionaries (one with the sequences filtered out 
                                     and another with the sequences kept in.)

    '''

    filtered_in = {}
    filtered_out = {}
    size_filter, char_filter, id_filter_out, id_filter_in = 0, 0, 0, 0
    for fasta_id in sequences:
        filter_activated = False
        
        seq = sequences[fasta_id].upper()
        seq_len = len(seq)
        
        # filter by sequence size
        if filter_by['seq_size']:
            if seq_len < min_seq_len:
                filtered_out[fasta_id] = sequences[fasta_id]
                size_filter += 1
                filter_activated = True

        # filter by sequence characters
        if filter_by['seq_char']:
            if all(char in seq for char in forbidden):
                filtered_out[fasta_id] = sequences[fasta_id]
                char_filter += 1
                filter_activated = True

        # filter by ID expressions (filter out or in)
        # filter out is prioritary
        if filter_by['id_list_out']:
            if any(s in fasta_id for s in id_filters):
                filtered_out[fasta_id] = sequences[fasta_id]
                id_filter_out += 1
                filter_activated = True
        elif filter_by['id_list_in']:
            if not any(s in fasta_id for s in id_filters):
                filtered_out[fasta_id] = sequences[fasta_id]
                id_filter_in += 1
                filter_activated = True
        
        # if no filters out are applied, apply filter in
        if not filter_activated:
            filtered_in[fasta_id] = sequences[fasta_id]
        filter_activated = False
            
    # Stats
    len_sequences = len(sequences)
    len_filtered = len(filtered_in)
    len_removed = len_sequences - len_filtered
    pct = len_filtered * 100.00 / len_sequences 
    
    # Prints
    if explain:
        print("\nFiltering results:")
        print("Initially, there were {} sequences.".format(len_sequences))
        print('Filtered out {} sequences with less than {} amino acids'.format(
            size_filter, min_seq_len))
        print('Filtered out {} sequences containing {} chars'.format(
            char_filter, ','.join(forbidden)))
        print('Filtered out {} sequences containing ids in list'.format(
            id_filter_out))
        print('Filtered out {} sequences not containing ids in list'.format(
            id_filter_in))
        print(">> Filtered out {} of {} sequences.".format(len_removed, 
                                                           len_sequences))
        print(">> Kept {} ({:.2f} %) sequences.".format(len_filtered, pct))
    
    if filtered_out:
        return {'in' : filtered_in, 'out': filtered_out}
    else:
        return filtered_in


def insert_newlines(string, every=64):
    '''
    Insert periodic newlines in a string.

    Parameters
    ----------
    string : STR
        String that will be modified.
    every : INT, optional
        Periodic position in which the string will be modified.
        The default is 64.

    Returns
    -------
    STR
        The modified string, with periodic newlines.

    '''
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))


def fasta_dict2file(fasta_dict, output_dir, output_file):
    '''
    Exports a fasta dictionay to a fasta file.

    Parameters
    ----------
    fasta_dict : DICT
        Dictionary containing fasta_ids as keys and sequences as values.
    output_dir : PATH
        Path of the output directory.
    output_file : STR
        Name of the output file.

    Returns
    -------
    None.

    '''
    
    output_path = output_dir / output_file
    print("\nExporting sequences to file {}".format(output_path))
    final_string = ''
    for fasta_id in fasta_dict:
        sequence = insert_newlines(fasta_dict[fasta_id], every=64)
        this_entry = ">{}\n{}\n".format(fasta_id, sequence)
        final_string += this_entry
    with open(output_path, 'w') as out:
        out.write(final_string)
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

# File paths
fasta_file = input_dir / 'uniprot-reviewed yes+AND+proteome up000005640.fasta'
list_file = input_dir / 'human_proteines_list.txt'

# Read files
fasta = fasta_parser(fasta_file, forbidden_lines=['',' '], verbose=True)
ids_list = list_parser(list_file, forbidden_lines=['',' '], verbose=True)

# Filter fasta
filtered_fastas = filter_sequences(fasta, 
                                  min_seq_len=35,
                                  forbidden=['B','J','O','U','X','Z'],
                                  id_filters=ids_list,
                                  filter_by={'seq_size':True,
                                             'seq_char':True,
                                             'id_list_in':True,
                                             'id_list_out':False},
                                  explain=True)

filtered_in  = filtered_fastas['in']
filtered_out = filtered_fastas['out']

# Export 
fasta_dict2file(filtered_in, output_dir=output_dir, 
                output_file='filtered_human.fasta')

