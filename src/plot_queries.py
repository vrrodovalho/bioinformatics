# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:48:34 2019

@author: vrrodovalho

This script contains a function for generating a line plot of number of 
publications in a year series from at least 2 csv files generated from
Pubmed queries.

"""

import os
import sys
import pathlib
import glob
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from textwrap import wrap


def plot_queries(query_dir, output_dir, file_name='plot.tiff', id_col='Year',
                 legend_pos='below', drop_values=[]):
    '''
    Generates a line plot from at least 2 csv files of Pubmed queries.

    Parameters
    ----------
    query_dir : STR
        The directory with the csv files from Pubmed queries.
    output_dir : STR
        The directory where the output graph will be generated.
    file_name : STR, optional
        The name of the output file. The default is 'plot.tiff'.
    id_col : STR, optional
        The name of the column containing ID values. The default is 'Year'.
    legend_pos : STR, optional
        The position of the legend. Possible values: ['below','right'].
        The default is 'below'.
    drop_values : LIST, optional
        List of ID values to be removed from the graph. The default is [].

    Returns
    -------
    new_merged : DataFrame
        A DataFrame containining all the data used to plot the graph.

    '''
    
    dfs = []
    # loop through csv files in query_dir and read dataframes
    for filename in glob.glob(os.path.join(query_dir, '*.csv')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            df = pd.read_csv(f, header=1)
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            header = f.readline().strip().split('Search query: ')[1]
        df = df.rename(columns={'Count': header})
        dfs.append(df)
    
    # merge all dataframes based on id_col
    idexed_dfs = [df.set_index(id_col) for df in dfs]
    merged = pd.concat(idexed_dfs, axis=1).reset_index() 
    
    # fill missing years
    min_year = merged[id_col].min()
    max_year = merged[id_col].max()
    new_df = pd.DataFrame({id_col:np.arange(min_year, max_year)})
    new_merged = new_df.merge(merged, on=id_col, how='outer')
    new_merged = new_merged.fillna(0).sort_values(by=id_col, ascending=True)
    new_merged = new_merged[~new_merged[id_col].isin(drop_values)]
    df = new_merged

    # prepare legend positions
    legend_pos_map = {'below' : {'dims' : (0, -0.55), 'wrap' : 80 }, 
                      'right' : {'dims' : (1.05, 0.5), 'wrap' : 30 }}
    bbox_to_anchor = legend_pos_map[legend_pos]['dims']
    legend_wrap = legend_pos_map[legend_pos]['wrap']

    # plot
    sns.set_style('darkgrid')
    random.seed(22)
    fig, ax = plt.subplots()
    labels = [ '\n'.join(wrap(l, legend_wrap)) \
              for l in df.columns if l != id_col ]
    df.plot.line(x='Year', ax=ax)
    ax.legend(labels, 
              title="Search query", 
              fontsize=8,
              title_fontsize=10,
              bbox_to_anchor=bbox_to_anchor, 
              loc='lower left', 
              borderaxespad=0.)
    ax.set_ylabel('Publications count')
    ax.set_xlabel('Year')
    
    # plt.subplots_adjust(left=0.1, right = 0.7)
    plt.show()
    
    fig.savefig(output_dir/file_name, 
                format='tiff', 
                dpi=300, 
                bbox_inches="tight")
    
    return new_merged

##############################################################################
    
# DIRECTORY SYSTEM
src_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
main_dir = os.path.dirname(src_dir)
root_dir = os.path.dirname(main_dir)
data_dir = pathlib.Path(main_dir) / 'data'
input_dir = pathlib.Path(data_dir) / 'input'
output_dir = pathlib.Path(data_dir) / 'output'
sys.path.insert(0, root_dir)

# DIRECTORY CONTAINING THE CSV FILES
query_dir = input_dir / 'queries'

# CALL PLOT FUNCTION
final_table = plot_queries(query_dir=query_dir, 
                           output_dir=output_dir, 
                           file_name='plot.tiff', 
                           id_col='Year',
                           legend_pos='below',
                           drop_values=[])

