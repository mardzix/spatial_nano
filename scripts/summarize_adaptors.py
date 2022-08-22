# Read structure:
#                                GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC -> Read starts here
# 5 AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC NNNNNNNN ATCCACGTGCTTGAGAGGCCAGAGCATTCG NNNNNNNN AGATGTGTATAAGAGACAGCATCGGCGTACGACT NNNNNN AGATGTGTATAAGAGACAG 3
#                                                                    * Spatial bcd 1 = bcd1           * Spatial barcode 2 = bcd3             * Modality barcode = bcd3

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import argparse
import regex
import yaml
import sys


ap = argparse.ArgumentParser()
ap.add_argument("-i",  "--input", required=True, help="input file")
ap.add_argument("-o",  "--output", help="Output file with per read adaptor statistics")
args = ap.parse_args()

spacer1     = 'ATCCACGTGCTTGAGAGGCCAGAGCATTCG'
spacer2     = 'AGATGTGTATAAGAGACAGCATCGGCGTACGACT'
spacer3     = 'AGATGTGTATAAGAGACAGCATCGGCGTACGACT.*AGATGTGTATAAGAGACAG'
spacer3_len = 34 + 8 + 19
spacer_P7   = 'CTGTCTCTTATACACATC'

statistics = {#'Spacer_P7 found': 0,
              'Spacer_1 not found' : 0,
              'Spacer_2 not found' : 0,
              'Spacer_3 not found' : 0,
              'Barcode wrong length': 0,
              'Reads with modality barcode': 0,
              'Read too short': 0,
              'Good read in modality': 0
              }

def find_seq(pattern, DNA_string, nmismatch=2):
    for n in range(0,nmismatch + 1):
        r = regex.compile('({0}){{e<={1}}}'.format(pattern, n))
        res = r.finditer(DNA_string)
        hit = [x.start() for x in res]
        if len(hit) == 0:
            continue
        if len(hit) > 1:
            return None
        if len(hit) == 1:
            return int(hit[0])
    return None

read_count = 0
with gzip.open(args.input, "rt") as in_handle_R1, open(args.output, "wt") as out_handle:
    for title, seq, qual in FastqGeneralIterator(in_handle_R1):
        read_count += 1
        # if read_count == 10000: # Debuging purposes
        #     break

        # Find spacers
        sp1   = find_seq(spacer1, seq)
        sp2   = find_seq(spacer2, seq)
        sp3   = find_seq(spacer3, seq)
        sp_P7 = find_seq(spacer_P7, seq,nmismatch=4)

        bcd1 = bcd2 = bcd3 = None

        if sp1:
            bcd1     = seq[:sp1]
            bcd1_q   = qual[:sp1]

        if sp1:
            bcd2     = seq[sp1 + len(spacer1):sp1 + len(spacer1) + 8]
            bcd2_q   = qual[sp1 + len(spacer1):sp1 + len(spacer1) +8]
        if sp2:
            bcd3     = seq[sp2 + len(spacer2):sp2 + len(spacer2) + 8]
            bcd3_q   = qual[sp2 + len(spacer2):sp2 + len(spacer2) + 8]

        if sp3 != None:
            new_read   = seq[sp3 + spacer3_len:sp_P7]    # If sp_P7 is not found then: sp_P7 = None
            new_read_q = qual[sp3 + spacer3_len:sp_P7]   # If sp_P7 is not found then: sp_P7 = None

        if sp3 != None:
            sp3_start = sp3 + 34 + 8
        else:
            sp3_start = None
        # sys.stdout.write(seq + '\n')
        # sys.stdout.write('read_{rc}\t{sp1}\t{sp2}\t{sp3}\t{sp_P7}\n'.format(rc=read_count, sp1=sp1, sp2=sp2, sp3=sp3_start,sp_P7=sp_P7))
        out_handle.write('read_{rc}\t{sp1}\t{sp2}\t{sp3}\t{sp_P7}\t{bcd1}\t{bcd2}\t{bcd3}\n'.format(rc = read_count,sp1 = sp1, \
                                                                                                    sp2 = sp2, sp3 = sp3_start, sp_P7 = sp_P7, \
                                                                                                    bcd1 = bcd1, bcd2 = bcd2, bcd3 = bcd3))


    yaml.dump(statistics, sys.stderr)
