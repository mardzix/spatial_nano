from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import gzip

ap = argparse.ArgumentParser()
ap.add_argument("--input_R1", required=True, help="input file after filtering R1")
ap.add_argument("--input_R2", required=True, help="input file for R2")
ap.add_argument("--output", required=True, help="output file read list")
args = ap.parse_args()

with gzip.open(args.input_R1, "rt") as in_handle1, \
     gzip.open(args.input_R2, 'rt') as in_handle2, \
     open(args.output, "wt") as out_handle:
    reads_R1 = []
    reads_R2 = []
    counter1 = 0
    counter2 = 0
    for title, seq, qual in FastqGeneralIterator(in_handle1):
        reads_R1.append(title.split(" ")[0])
        # counter1 += 1
        # if counter1 == 100000:
        #     break
    for title, seq, qual in FastqGeneralIterator(in_handle2):
        reads_R2.append(title.split(" ")[0])
        # counter2 += 1
        # if counter2 == 100000:
        #     break
    for read in list(set(reads_R1 + reads_R2)):
        if read in reads_R1 and read in reads_R2:
            out_handle.write(read + '\n')




