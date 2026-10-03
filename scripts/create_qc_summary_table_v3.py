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

vars.add_argument('--fasta',dest='assembly_fasta', type=str, required=False,default='', help='assembly fasta file')

vars.add_argument('--gfastats',dest='gfastats_file', type=str, required=False,default='', help='stats.txt file produce by gfastats')

vars.add_argument('--compleasm',dest='compleasm_file', type=str, required=False,default='', help='stats.txt file produced by compleasm')

vars.add_argument('-o','--output_prefix',dest='output_prefix', type=str, required=False,default='', help='prefix for the output file')

vars.add_argument('-od','--out_dir',dest='output_directory', type=str, required=False,default=None, help='name of directory for output. will be stored in current working directory if none given')

vars.add_argument('--assem_name',dest='assembly_name', type=str, required=False,default='SpeciesX_Haplomex', help='Name for the assembly i.e Haplome1')

vars.add_argument('-c',dest='chrs', type=int, required=False,default=8, help='Expected number of chromsomes in the assembly. Default = 8')

args = vars.parse_args()

# by default directory to write output to is current wd
outdir = os.getcwd()

if args.output_directory:
    outdir = args.output_directory # set output dir

# check if dir given exists, if not create one
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# configure the logger we will be using
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.FileHandler(f"{outdir}/create_qc_summary_tables.log",mode='w')],datefmt='%H:%M:%S')

# create a logger
Mainlogger = logging.getLogger(__name__)

# set main summary table
MainSummaryTable = {'':[],
                    args.assembly_name:[]}

with pd.ExcelWriter(f'{outdir}/{args.output_prefix}_MasterSpreadsheet.xlsx', engine='openpyxl') as ExcelOut:

    # if fasta file provided
    if args.assembly_fasta:

        with open(args.assembly_fasta,'r') as AssemblyFastaFile:

            # create dict for assembly info
            AssemblyInfo = {'ChromID':[],
                            'ChromSize':[]}


            # store the fasta data for the assembly in a list but do not remove new lines
            AssemblyFastaData = AssemblyFastaFile.readlines()

            # chromosomes and their sequence will be stored here
            AssemblyFastaRefDict = fasta_to_dict(AssemblyFastaData)

            # get the length of all chromosomes in the assembly
            AssemblyStatistics = {}

            # sort dict based on size of chromosomes
            AssemblyFastaRefDict = dict(sorted(AssemblyFastaRefDict.items(),key=lambda name_size:name_size[1].sequence_length,reverse=True))

            for Chr in AssemblyFastaRefDict:
                Mainlogger.info(f'Chromosome {Chr.rstrip()} has length {AssemblyFastaRefDict[Chr].sequence_length}')

            AssemblyInfo['ChromID'].extend(list(AssemblyFastaRefDict.keys())[0:args.chrs])
            AssemblyInfo['ChromSize'].extend([f"{record.sequence_length:,}" for record in list(AssemblyFastaRefDict.values())[0:args.chrs]])

            AssemblyInfoDf = pd.DataFrame(AssemblyInfo)

            AssemblyInfoDf.to_excel(excel_writer=ExcelOut,sheet_name='AssemblyInfo',index=False)
            AssemblyInfoDf.to_csv(path_or_buf=f'{outdir}/{args.output_prefix}_AssemblyInfo.tsv',sep='\t',index=False)


    # if gfastats file provided
    if args.gfastats_file:

        with open(args.gfastats_file,'r') as GfastatsFile:

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

            GfastatsDict = dict(zip(GfastatsInfo[0::2],GfastatsInfo[1::2]))

            # uncomment to know exact indeces of every line of data in gfastats output file
            # for i,value in enumerate(GfastatsData):
            #     print(f'Index {i} = {value}')


            # create dict for main summary table    
            MainSummaryTable[''].extend(['Number of Scaffolds','Total Length (bp)','Scaffold N50 (bp)','Largest Scaffold (bp)','Number of Contigs','Contig N50','Largest Contig (bp)','Gaps','GC Content(%)'])
            
            MainSummaryTable[args.assembly_name].append(GfastatsDict['scaffolds'])
            MainSummaryTable[args.assembly_name].append(GfastatsDict['Total scaffold length'])
            MainSummaryTable[args.assembly_name].append(GfastatsDict['Scaffold N50'])
            MainSummaryTable[args.assembly_name].append(GfastatsDict['Largest scaffold'])
            MainSummaryTable[args.assembly_name].append(GfastatsDict['contigs'])
            MainSummaryTable[args.assembly_name].append(GfastatsDict['Contig N50'])
            MainSummaryTable[args.assembly_name].append(GfastatsDict['Largest contig'])
            MainSummaryTable[args.assembly_name].append(GfastatsDict['gaps'])
            MainSummaryTable[args.assembly_name].append(GfastatsDict['GC content %'])
    
    # if compleasm data provided
    if args.compleasm_file:

        with open(args.compleasm_file) as CompleasmFile:

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
            MainSummaryTable[''].append('BUSCOs(%)')
            MainSummaryTable[args.assembly_name].append(BuscoData)



    # summary table to df
    MainSummaryTableDf = pd.DataFrame(MainSummaryTable)

    MainSummaryTableDf.to_excel(excel_writer=ExcelOut,sheet_name='Main',index=False)
    MainSummaryTableDf.to_csv(path_or_buf=f'{outdir}/{args.output_prefix}_MasterSpreadsheet.tsv',sep='\t',index=False)
