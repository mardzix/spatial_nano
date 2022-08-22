# Wrapper around the nano-CUT&Tag spatial data analysis pipeline
# 07/07/2022
# The pipeline takes fastq files with multiplexed and spatial nano-CUT&Tag data as input
# It performs:
#   1. Demultiplexing of the fastq files
#   2. Run cellranger-atac
#   3. Construct seurat object with pixel x peak matrix

import argparse
import glob
import os, sys
import yaml
import fnmatch
import itertools

parser = argparse.ArgumentParser(description='Preprocessing pipeline for single-cell nano-CUT&Tag data analysis')
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fastq',  type=str, required=True, help='Path to folder with fastq files')
parser.add_argument('--sample',  type=str, required=False,default='Audodetect', help='Sample name override')
parser.add_argument('--barcodes',  type=str, required=True,nargs="+", help='list of space separated barcodes to be used for demultiplexing [e.g. ATAGAGGC TATAGCCT CCTATCCT]')
parser.add_argument('--modalities',  type=str, required=True,nargs="+", help='list of space separated modalities corresponding to the barcodes [e.g. ATAC H3K27ac H3K27me3]')
parser.add_argument('--cellranger_ref', type=str,default = '/data/ref/cellranger-atac/refdata-cellranger-atac-mm10-2020-A-2.0.0/', help = 'path to cellranger reference folder')
parser.add_argument('--threads',  type=int, default = 1, help='Number of threads')
parser.add_argument('--genome',  type=str, default = 'mm10', help='Genome to use for Gene activity scores [only mm10 supported for now]')
parser.add_argument('--tempdir',  type=str, default = './temp', help='Path to temp directory for sort command')
parser.add_argument('--condaenv',  type=str, default = False, help='Conda environment to use for analysis')
parser.add_argument('--snakeargs', default = " ",  type=str, nargs="+", help='Optional arguments to be passed to snakemake, must be formated --snakeargs="--PLACE_ARGS_HERE --MORE_ARGS" [e.g. --snakeargs="--dryrun --printshellcmds"')
args = parser.parse_args()

fastq_ext = '.fq.gz'

def get_sample_id_from_fastq_folder(folder):
    folder      = os.path.abspath(folder)
    sample_name = folder.split('/')
    sample_name = sample_name[len(sample_name)-1]
    return(sample_name)

def match_barcode_to_modality(modalities, barcodes):
    if len(modalities) != len(barcodes):
        sys.stderr.write("*** ERROR: Number of modalities does not match the numbe of barcodes\n")
        sys.stderr.write("Modalities: {}\nBarcodes:{}\n".format(modalities,barcodes))
        sys.exit()
    return {modalities[i]: barcodes[i] for i,x in enumerate(modalities)}

#####################################
# Detect sample name
if args.sample == 'Audodetect':
    sample_id = get_sample_id_from_fastq_folder(args.fastq)
else:
    sample_id = args.sample
sys.stderr.write('\n*** Log: Sample ID: {} \n'.format(sample_id))
######################################

# Find fastq files within the args.fastq folder
files = [os.path.abspath(x) for x in glob.glob(args.fastq + "/*")]
patterns_R1 = [a+b for a in ['*_R1_*','*_1_*'] for b in ['.fq.gz','.fastq.gz'] for s in [sample_id]]
patterns_R2 = [a+b for a in ['*_R2_*','*_2_*'] for b in ['.fq.gz','.fastq.gz'] for s in [sample_id]]
files = {
            'R1': list(itertools.chain(*[fnmatch.filter(files,p) for p in patterns_R1])) ,
            'R2': list(itertools.chain(*[fnmatch.filter(files,p) for p in patterns_R2]))
}
sys.stderr.write('\n*** Log: Files found in the folder: \n    R1: {}\n    R2: {}\n'.format(files['R1'],files['R2']))


# Create config
SNAKEFILE = os.path.dirname(os.path.realpath(__file__)) + "/workflow/Snakefile"
modality_barcode_dict = match_barcode_to_modality(args.modalities, args.barcodes)

CONFIG = {}
CONFIG['cellranger_ref'] = args.cellranger_ref
CONFIG['input_folder']   = args.fastq
CONFIG['fastq']          = files
CONFIG['sample']         = sample_id
CONFIG['modalities']     = modality_barcode_dict
CONFIG['genome']         = args.genome
CONFIG['tempdir']        = args.tempdir
CONFIG['condaenv']       = args.condaenv

with open(sample_id + '_config.yaml', 'w') as outfile:
    sys.stderr.write("\n*** Log: Config dump: \n")
    yaml.dump(CONFIG,sys.stderr)
    yaml.dump(CONFIG, outfile, default_flow_style=False)

pipeline_cmd = "snakemake --snakefile {snakefile} --cores {threads} --configfile {configfile} -p {pass_args}".format(
    snakefile=SNAKEFILE,
    threads = args.threads,
    configfile = sample_id + '_config.yaml',
    pass_args = " ".join(args.snakeargs))

sys.stderr.write("\n*** Log: Pipeline command:\n {}\n ".format(pipeline_cmd))
os.system(pipeline_cmd)