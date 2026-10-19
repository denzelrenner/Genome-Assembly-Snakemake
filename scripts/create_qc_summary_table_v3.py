#!/usr/bin/env python3

import os
import argparse
import logging
import pandas as pd
import json
import sys
from PowerhouseOfTheCell import *

# This script will take in the output files from gfastats and busco
# Input
# gfastats stats.txt file
# busco summary json file
# assembly fasta file

# Output
# excel spreadsheet with the data parsed


vars = argparse.ArgumentParser(description='Create summary statistics from gfastats and busco')

vars.add_argument('-s','--samples',dest='sample',nargs='+', type=str, required=True, help='sample(s) to process stats for')

vars.add_argument('-nhap',dest='haplomes',nargs='+', type=int, required=True, help='number of haplomes after running hifiasm')

vars.add_argument('-o','--output_prefix',dest='outputPrefix', type=str, required=False,default='output', help='prefix for the output file')

vars.add_argument('-od','--out_dir',dest='output_directory', type=str, required=True,default=None, help='name of directory for output. will be stored in current working directory if none given')

vars.add_argument('-c',dest='chrs', type=int, required=False,default=8, help='Expected number of chromsomes in the assembly. Default = 8')

args = vars.parse_args()

# by default directory to write output to is current wd
outdir = os.getcwd()

if args.output_directory:
    outdir = args.output_directory # set output dir

# check if dir given exists, if not create one
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# define functions
def get_chromosome_sizes(haplome_name:str,file:str):

    # create dict for assembly info
    AssemblyInfo = {haplome_name:[]}

    # chromosomes and their sequence will be stored here
    AssemblyFastaRefDict = load_fasta(file)

    # sort dict based on size of chromosomes
    AssemblyFastaRefDict = dict(sorted(AssemblyFastaRefDict.items(),key=lambda name_size:name_size[1].sequence_length,reverse=True))

    # add chrom id and its size separated by colon
    AssemblyInfo[haplome_name].extend([f"{record.header} ---> {record.sequence_length:,}" for record in list(AssemblyFastaRefDict.values())])

    return pd.DataFrame(AssemblyInfo)


def gfastats_to_df(haplome_name:str,file:str):

    with open(file,'r') as GfastatsFile:

        # set main summary table
        MainSummaryTable = {'':['Number of Scaffolds','Total Length (bp)','Scaffold N50 (bp)','Largest Scaffold (bp)','Number of Contigs','Contig N50','Largest Contig (bp)','Gaps','GC Content(%)'],haplome_name:[]}

        ### Parse gfastats data
        GfastatsData = GfastatsFile.readlines()

        # remove index 25 with base composition and remove index 0 becausse they dont contain important information
        del GfastatsData[25]
        del GfastatsData[0]

        GfastatsInfo = []
        for Line in GfastatsData:

            # remove bad characters
            Info = Line.split(':')[0].lstrip('#').lstrip().rstrip()
            UntypedNum = Line.split(':')[1].lstrip().rstrip()
            Data =  f"{float(UntypedNum):,}" if '.' in UntypedNum else f"{int(UntypedNum):,}"
            GfastatsInfo.append(Info)
            GfastatsInfo.append(Data)

        # create dict from parsed data
        GfastatsDict = dict(zip(GfastatsInfo[0::2],GfastatsInfo[1::2]))

        # update dict for main summary table    
        MainSummaryTable[haplome_name].append(GfastatsDict['scaffolds'])
        MainSummaryTable[haplome_name].append(GfastatsDict['Total scaffold length'])
        MainSummaryTable[haplome_name].append(GfastatsDict['Scaffold N50'])
        MainSummaryTable[haplome_name].append(GfastatsDict['Largest scaffold'])
        MainSummaryTable[haplome_name].append(GfastatsDict['contigs'])
        MainSummaryTable[haplome_name].append(GfastatsDict['Contig N50'])
        MainSummaryTable[haplome_name].append(GfastatsDict['Largest contig'])
        MainSummaryTable[haplome_name].append(GfastatsDict['gaps'])
        MainSummaryTable[haplome_name].append(GfastatsDict['GC content %'])

        return pd.DataFrame(MainSummaryTable)

def compleasm_to_df(haplome_name:str,file:str):

    with open(file) as CompleasmFile:

        # set main summary table
        MainSummaryTable = {'':['BUSCOs(%)'],haplome_name:[]}

        ### Parse Compleasm data
        # store busco data
        CompleasmData = CompleasmFile.readlines()[1:]
        TotalBuscoGenes = int(CompleasmData[-1].split(':')[-1])
        CompleasmData = CompleasmData[:-1] # remove last line
        CompleasmDict = {} 
        for Line in CompleasmData:
            Frequency = Line.split(',')[-1].lstrip().rstrip()
            Classification,Percentage = Line.split(',')[0].split(':')
            CompleasmDict[Classification] = (Percentage,Frequency)
        
        # add complete buscos to the dict
        NoCompleteBuscos = int(CompleasmDict['S'][1])+int(CompleasmDict['D'][1])
        # CompleasmDict['C'] = (f"{round(NoCompleteBuscos/TotalBuscoGenes,4)*100}%",f"{NoCompleteBuscos}")
        CompleteBuscoPercentage = float(CompleasmDict['S'][0].rstrip('%')) + float(CompleasmDict['D'][0].rstrip('%'))
        CompleasmDict['C'] = (f"{CompleteBuscoPercentage}%",f"{NoCompleteBuscos}")
        BuscoData = f'''C:{CompleasmDict['C'][0]}; {CompleasmDict['C'][1]}\nS:{CompleasmDict['S'][0]}; {CompleasmDict['S'][1]}\nD:{CompleasmDict['D'][0]}; {CompleasmDict['D'][1]}\nF:{CompleasmDict['F'][0]}; {CompleasmDict['F'][1]}\nI:{CompleasmDict['I'][0]}; {CompleasmDict['I'][1]}\nM:{CompleasmDict['M'][0]}; {CompleasmDict['M'][1]}\nN:{TotalBuscoGenes}'''

        # add info to main summary table
        MainSummaryTable[haplome_name].append(BuscoData)

    return pd.DataFrame(MainSummaryTable)


# run different analysis depending on if we have a single sample or many
if len(args.sample) == 1:

    # merged qc stats and chrom sizes
    toMergeQcStats = []
    toMergeChromSizes = []

    for i in range(1,args.haplomes[0]+1):

        # get haplome name
        haplome = f'{args.sample[0]}.hap{i}'

        # get chrom sizes df
        chromSizes = get_chromosome_sizes(haplome_name=haplome,file=f'{args.sample[0]}/s4_run_hifiasm/{args.sample[0]}.hap{i}.fa')

        # get gfastats df
        gfastats = gfastats_to_df(haplome_name=haplome,file=f'{args.sample[0]}/s4_run_hifiasm/quality_control/{args.sample[0]}.hap{i}/gfastats/stats.txt')

        # get compleams df
        compleasm = compleasm_to_df(haplome_name=haplome,file=f'{args.sample[0]}/s4_run_hifiasm/quality_control/{args.sample[0]}.hap{i}/compleasm/summary.txt')
       
        # combine stats
        qcStats = concat_multiple_dataframes([gfastats,compleasm])

        # add to dfs to merge 
        toMergeQcStats.append(qcStats)
        toMergeChromSizes.append(chromSizes)

    # merge dfs and then output
    mergedQcStats = merge_multiple_dataframes(dfs=toMergeQcStats,on='')
    mergedChromSizes = merge_multiple_dataframes(dfs=toMergeChromSizes,how='outer',left_index=True,right_index=True)

    # create df sheet name output
    df_sheetname = [(mergedQcStats,'qc stats'),(mergedChromSizes,'scaffold info')]
    
    # output 
    many_df_to_excel(df_sheetName=df_sheetname,outputDir=outdir,outputFilePrefix=args.outputPrefix,col_width=40)

elif len(args.sample) > 1:

    # merged qc stats and chrom sizes
    toMergeQcStats = []
    toMergeChromSizes = []

    # make pairs of sample and haplome
    sample_haplome = list(zip(args.sample,args.haplomes))
    
    for s,h in sample_haplome:

        for i in range(1,h+1):

            # get haplome name
            haplome = f'{s}.hap{i}'

            # get chrom sizes df
            chromSizes = get_chromosome_sizes(haplome_name=haplome,file=f'{s}/s4_run_hifiasm/{s}.hap{i}.fa')

            # get gfastats df
            gfastats = gfastats_to_df(haplome_name=haplome,file=f'{s}/s4_run_hifiasm/quality_control/{s}.hap{i}/gfastats/stats.txt')

            # get compleams df
            compleasm = compleasm_to_df(haplome_name=haplome,file=f'{s}/s4_run_hifiasm/quality_control/{s}.hap{i}/compleasm/summary.txt')
        
            # combine stats
            qcStats = concat_multiple_dataframes([gfastats,compleasm])

            # add to dfs to merge 
            toMergeQcStats.append(qcStats)
            toMergeChromSizes.append(chromSizes)

    # merge dfs and then output
    mergedQcStats = merge_multiple_dataframes(dfs=toMergeQcStats,on='')
    mergedChromSizes = merge_multiple_dataframes(dfs=toMergeChromSizes,how='outer',left_index=True,right_index=True)

    # create df sheet name output
    df_sheetname = [(mergedQcStats,'qc stats'),(mergedChromSizes,'scaffold info')]
    
    # output 
    many_df_to_excel(df_sheetName=df_sheetname,outputDir=outdir,outputFilePrefix=args.outputPrefix,col_width=40)


