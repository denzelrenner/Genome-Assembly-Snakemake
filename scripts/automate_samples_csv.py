import argparse
import os
import pandas as pd

vars = argparse.ArgumentParser(description='Create template samples.csv file given a directory of Hifi reads')

vars.add_argument('-i','--hifi',dest='hifiDir', type=str, required=True, help='absolute path to directory with hifi reads')

vars.add_argument('-od','--out_dir',dest='output_directory', type=str, required=False,default=None, help='name of directory for output. will be stored in current working directory if none given')

vars.add_argument('--exclude',dest='samplesToExclude',nargs='+',type=str, required=False,default=None, help='samples to exclude from the input directory, can be space separated list e.g sampleA sampleB sampleC; or can be path to txt file e.g ./exclude.txt')

vars.add_argument('--csv-columns',dest='csvColumns',nargs='+', type=str, required=False,default=['useHiC','False','ploidy','2','chrom','8','haplomeSize','200'], help='column names and default values for output csv. Do not include default column for samples e.g useHiC False ploidy 2')

args = vars.parse_args()

# by default directory to write output to is current wd
outdir = os.getcwd()

if args.output_directory:
    outdir = args.output_directory # set output dir

# check if dir given exists, if not create one
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# load in hifi reads dir and get sample name by removing .hifi_reads
samples = [i.split('.hifi_reads')[0] for i in os.listdir(args.hifiDir)]

# check for sample to exclude
if args.samplesToExclude:

    # check if file
    if len(args.samplesToExclude) == 1 and args.samplesToExclude[0].endswith('.txt'):

        # load in txt 
        excludedSamples = pd.read_csv(args.samplesToExclude[0],header=None,names=['sample'])['sample'].to_list()
        
        # subtract from samples
        samples = list(set(samples) - set(excludedSamples))

    # otherwise it should be space separated list so just subtract from samples
    else:
        samples = list(set(samples) - set(args.samplesToExclude))


# create df as a dictionary, merge dict of only samples, and other column values;  use zip to create default
dfAsDict = {**{'sample':samples},**{colname:[defaultValue]*len(samples) for colname,defaultValue in zip(args.csvColumns[0::2],args.csvColumns[1::2])}}

# create output df
df = pd.DataFrame(data=dfAsDict)

# output csv
df.to_csv(path_or_buf=os.path.join(outdir,'samples.csv'),index=False,sep=',')


