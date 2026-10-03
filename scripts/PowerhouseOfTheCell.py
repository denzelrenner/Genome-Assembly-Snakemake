#!/usr/bin/env python3
import subprocess
import os
from functools import reduce
import numpy as np
import pandas as pd

class FastaRecord():

    def __init__(self,header:str,array_of_sequences:str):

        self.header = header.rstrip('\n')

        # the sequence is currently a list of sequences so we need to concatenate them together
        self.sequence = ''.join([Line.rstrip() for Line in array_of_sequences])

        # get the sequence length
        self.sequence_length = len(self.sequence)

        # when outputting the sequence we want to preserve the original structure of the input file including new lines
        self.output_sequence = ''.join(array_of_sequences)
        self.output_header = ">" + header.rstrip('\n') + '\n'

        # get the number of each base
        self.A_count,self.C_count,self.T_count,self.G_count = find_ACTG_content(sequence=self.sequence)
        self.A_percentage = (self.A_count/self.sequence_length) * 100
        self.C_percentage = (self.C_count/self.sequence_length) * 100
        self.T_percentage = (self.T_count/self.sequence_length) * 100
        self.G_percentage = (self.G_count/self.sequence_length) * 100

    # function to output record, i.e write header and sequence to a single file object, file must be opened prior to using this method
    def write_record_to_file(self,file_object) -> None:

        # write output
        file_object.write(self.output_header)
        file_object.writelines(self.output_sequence)


def create_multiline_fasta(sequence:str,characters_per_line:int=66) -> list:

    output = []

    for i in range(0,len(sequence),characters_per_line):

        output.append(f"{sequence[i:i+characters_per_line]}\n")
    
    return output

def subset_fasta(sequence:str,start:int,end:int) -> str:

    return sequence[start:end]
    
# useful functions
def find_ACTG_content(sequence:str) -> tuple:

    return sequence.count('A'),sequence.count('C'),sequence.count('T'),sequence.count('G')

# create dict with header name as key and Fasta Record object as value
def fasta_to_dict(FastaData:list) -> dict:

    # scaffold and their sequence will be stored here
    RefDict = {}

    # current seq is a list because we will be storing sequences without removing the new line, so we want each header to have a 'group' of sequences
    CurrentSeq = []
    CurrentHeader = FastaData[0]

    # generate dictionary with gene names and fasta files
    for i in range(1,len(FastaData)):

        # find what the next line is
        Line = FastaData[i]

        # find if the line begins with a gene
        if Line.startswith('>'): 

            # if line begins with a gene and there is some sequence for the previous gene that was found, add it to the dictionary reset the sequence builder
            if CurrentSeq:
                RefDict[CurrentHeader.replace('>','').rstrip()]=FastaRecord(header=CurrentHeader.replace('>','').rstrip(),array_of_sequences=CurrentSeq)
                CurrentHeader=Line

            CurrentSeq = []

        # if a line starts with the last line/sequence in the file, add it to the current sequence builder and then add the current gene to the ref dictionary
        elif i == len(FastaData) - 1:
            CurrentSeq.append(Line)
            RefDict[CurrentHeader.replace('>','').rstrip()]=FastaRecord(header=CurrentHeader.replace('>','').rstrip(),array_of_sequences=CurrentSeq)

        # if none of the above conditions are met then keep building the sequence
        else:
            CurrentSeq.append(Line)

    return RefDict

# create dict with scaffold name as key
def fasta_to_dict_original(FastaData:list):

    # scaffold and their sequence will be stored here
    RefDict = {}

    # current seq is a list because we will be storing sequences without removing the new line, so we want each header to have a 'group' of sequences
    CurrentSeq = []
    CurrentHeader = FastaData[0]

    # generate dictionary with gene names and fasta files
    for i in range(1,len(FastaData)):

        # find what the next line is
        Line = FastaData[i]

        # find if the line begins with a gene
        if Line.startswith('>'): 

            # if line begins with a gene and there is some sequence for the previous gene that was found, add it to the dictionary reset the sequence builder
            if CurrentSeq:
                RefDict[CurrentHeader.replace('>','')]=CurrentSeq
                CurrentHeader=Line

            CurrentSeq = []

        # if a line starts with the last line/sequence in the file, add it to the current sequence builder and then add the current gene to the ref dictionary
        elif i == len(FastaData) - 1:
            CurrentSeq.append(Line)
            RefDict[CurrentHeader.replace('>','')]=CurrentSeq

        # if none of the above conditions are met then keep building the sequence
        else:
            CurrentSeq.append(Line)

    return RefDict

# function to load fasta file, return dictionary of fasta header and fasta recod
def load_fasta(fastaFile:str) -> dict:

    AssemblyFastaRefDict = None

    with open(fastaFile,'r') as assemblyfastafile:

        # store the fasta data for the assembly in a list but do not remove new lines
        AssemblyFastaData = assemblyfastafile.readlines()

        # chromosomes and their sequence will be stored here
        AssemblyFastaRefDict = fasta_to_dict(AssemblyFastaData)

    return AssemblyFastaRefDict

# function to create df of fasta record and length of entry
def fasta_length_to_df(fastaFile:str):

    refDict = load_fasta(fastaFile)

    return pd.DataFrame(data={'gene':[v.header for k,v in refDict.items()],'fastalength':[v.sequence_length for k,v in refDict.items()]})





class Gff3Record():

    def __init__(self,gff_line:str,info_col_sep:str=';'):

        # breakdown the line from the gff
        self.complete_line = gff_line.split('\t')

        self.chr = self.complete_line[0]
        self.source = self.complete_line[1]
        self.feature = self.complete_line[2]
        self.start = self.complete_line[3]
        self.end = self.complete_line[4]
        self.score = self.complete_line[5]
        self.strand = self.complete_line[6]
        self.phase = self.complete_line[7]
        self.info = self.complete_line[8]

        if ';' in self.info:
            self.info_parsed = {i.split('=')[0].rstrip():i.split('=')[1].rstrip() for i in self.info.split(info_col_sep) if '=' in i} # need to write something to handle if the info field has no =
        
        else:
            self.info_parsed = {'ID':self.info.rstrip()}

    # rebuild and return a tab separated complete line, this way if we change the value in a column it can be reflected correctly
    def rebuild_and_return_line(self):
        
        return '\t'.join([self.chr,self.source,self.feature,self.start,self.end,self.score,self.strand,self.phase,self.info])
    
    # return the full line as was given in the input data
    def return_line(self):
        
        return '\t'.join(self.complete_line)
    
    # return the full line with only a specific tag from the info column
    def return_line_tag_only(self,tag:str):

        info = self.info_parsed[tag]

        line = self.complete_line
        line[8] = info

        return '\t'.join(line)

    # print line as a bed file
    def return_line_as_bed(self):

        return '\t'.join([self.chr,self.start,self.end])
    
    # print line as bed file with info 

    def return_line_as_bed_with_tag(self,tag:str):

        info = self.info_parsed[tag]

        return '\t'.join([self.chr,self.start,self.end,info])
    
    def return_line_as_bed_no_tag(self):

        return '\t'.join([self.chr,self.start,self.end,self.info])


#### Terminal ####

# general useful functions
def terminal_pipe(cmd): 
    stdout,stderr = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()

    return (stdout.decode("utf-8").strip(' \n'),stderr.decode("utf-8").strip(' \n'))

# shell is True, so you can use commands like pipes etc
def terminal_pipe_st(cmd):
    
    stdout,stderr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()

    return (stdout.decode("utf-8").strip(' \n'),stderr.decode("utf-8").strip(' \n'))



#### Blobtools ####

# function to take index values in blobdir json
def index_blobdir_json(data:dict,key:str='values') -> dict:

    # create out dict
    out_dict = {}

    # store list
    vals = data[key]

    # index list
    for i in range(len(vals)):

        out_dict[i] = vals[i]
    
    return out_dict


#### Fasta Files ####

# function to subset an entire fasta file based on a list of records we want to keep
def subset_whole_fasta(FastaRefDict:dict,records:list) -> dict:

    filtered_fasta = {}

    for header,record in FastaRefDict.items():

        if record.header in records:
            filtered_fasta[record.header] = record

    return filtered_fasta

# output fasta files
def output_fasta_file(FastaRefDict:dict,output_directory:str,file_name:str) -> None:

    # make sure outdir exists
    os.makedirs(output_directory,exist_ok=True)

    # create outfile
    with open(f'{output_directory}/{file_name}','w') as outfile:

        # only one sequence in the dict but we can safely loop
        for header,record in FastaRefDict.items():

            # write output
            outfile.write(record.output_header)
            outfile.writelines(record.output_sequence)


# output fasta file as upper
def output_fasta_file_upper(FastaRefDict:dict,file_name:str,output_directory:str=os.getcwd(),reformat=False) -> None:

    # create outfile
    with open(f'{output_directory}/{file_name}','w') as outfile:

        # only one sequence in the dict but we can safely loop
        for header,record in FastaRefDict.items():

            # write output
            outfile.write(record.output_header)

            if reformat:
                outfile.writelines(create_multiline_fasta(record.sequence.upper()))

            elif not reformat:
                outfile.writelines(record.output_sequence.upper())





# pandas toools
# def dict_to_excel(myDict:dict,outPath:str):

#     MainSummaryTableDf = pd.DataFrame(myDict)

#     MainSummaryTableDf.to_excel(excel_writer=ExcelOut,sheet_name='Main',index=False)
#     MainSummaryTableDf.to_csv(path_or_buf=f'{outdir}/{args.output_prefix}_MasterSpreadsheet.tsv',sep='\t',index=False)


### pandas tools

# function to merge multiple dataframes, input is list of dataframe 
def merge_multiple_dataframes(dfs:list,on:str='',how:str='inner'):
    return reduce(lambda left,right: pd.merge(left,right,on=on,how=how),dfs)

# write df to excel and tsv, only single sheet in excel file
def write_to_tsv_and_excel(df:pd.DataFrame,outputDir:str,outputFilePrefix:str,excel=True,csv=True,seperator='\t'):

    # make dirs if it does not exist
    os.makedirs(outputDir,exist_ok=True)

    # write to tsv and excel file
    if csv:
        df.to_csv(f'{outputDir}/{outputFilePrefix}.tsv',sep=seperator,index=False)

    if excel:
        df.to_excel(f'{outputDir}/{outputFilePrefix}.xlsx',index=False)
    
## apply functions

# function to join columns, (for a given row) we take 4 columns for example and merge the values in these columns into a single one
def join_columns(row,join_string:str=''):
    return join_string.join(row)

