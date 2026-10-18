#!/usr/bin/env python3

import os
import argparse
import logging
import pandas as pd
from PowerhouseOfTheCell import *

# Input:
# a tsv from seqkit
# an expected haplome size
# an outdir and out prefix

# Output:
# excel spreadsheet of information

vars = argparse.ArgumentParser(description='Create summary statistics from seqkit analysis for an individual sample')

vars.add_argument('-s','--sample',dest='sample', type=str, required=True, help='sample to process stats for')

vars.add_argument('-o','--output_prefix',dest='outputPrefix', type=str, required=False,default='output', help='prefix for the output file')

vars.add_argument('-od','--out_dir',dest='output_directory', type=str, required=True,default=None, help='name of directory for output. will be stored in current working directory if none given')

vars.add_argument('-g','--genome_size',dest='genomeSize', type=int, required=True,default=210000000, help='expected haplome size (in bp)')

args = vars.parse_args()

# function to create data frame from seqkit tsv input
def create_seqkit_df(inputFilePath:str):

    # open tsv file
    df = pd.read_csv(inputFilePath,sep='\t',header=0)

    # add column for sequencing depth
    df['sequencing_depth'] = round(df['sum_len']/args.genomeSize,2)

    # add column for total sequenced dna in gigabases
    df['sum_len (in Gigabases)'] = df['sum_len']//(10**9)

    # add column for file size, divide number of bytes by how much is in one gigabit (GiB), and also gigabyte (Gb)
    df['file_size (in Gib)'] = round(df['file'].apply(os.path.getsize)/(1024**3),2)
    df['file_size (in Gb)'] = round(df['file'].apply(os.path.getsize)/(10**9),2)

    # add , to everything that is not a string, i.e the first 3 columns of output
    numcols = list(set(df.columns.values) - set(['file','format','type']))
    df[numcols] = df[numcols].map(lambda x:f'{x:,}')

    return df

# by default directory to write output to is current wd
outdir = os.getcwd()

if args.output_directory:
    outdir = args.output_directory # set output dir

# check if dir given exists, if not create one
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# first check if we have hifi, hic1, and hic2 files
hifi = (create_seqkit_df(f'{args.sample}/s2_input_data_qc/Hifi/{args.sample}.hifi.tsv'),'hifi') if os.path.exists(f'{args.sample}/s2_input_data_qc/Hifi/{args.sample}.hifi.tsv') else None
hic1 = (create_seqkit_df(f'{args.sample}/s2_input_data_qc/HiC1/{args.sample}.hic1.tsv'),'hic1') if os.path.exists(f'{args.sample}/s2_input_data_qc/HiC1/{args.sample}.hic1.tsv') else None
hic2 = (create_seqkit_df(f'{args.sample}/s2_input_data_qc/HiC2/{args.sample}.hic2.tsv'),'hic2') if os.path.exists(f'{args.sample}/s2_input_data_qc/HiC2/{args.sample}.hic2.tsv') else None

# array to store tuples of df and sheet name
dataframes_sheetnames = [i for i in [hifi,hic1,hic2] if i]

many_df_to_excel(dataframes_sheetnames,outputDir=outdir,outputFilePrefix=args.outputPrefix)
