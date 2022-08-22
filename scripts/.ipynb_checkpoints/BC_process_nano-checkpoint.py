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
import Levenshtein


ap = argparse.ArgumentParser()
ap.add_argument("-i",  "--input", required=True, help="input file")
ap.add_argument("-o1", "--output_R1", required=True, help="output file R1")
ap.add_argument("-o2", "--output_R2", required=True, help="output file R2")
ap.add_argument("-l",  "--log", default = 'log.yaml',  help="Log of output statistics")
ap.add_argument('-b',  "--barcode", required=True, help="")
args = ap.parse_args()

input_file_R1 = args.input
output_file_R1 = args.output_R1
output_file_R2 = args.output_R2

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
with gzip.open(input_file_R1, "rt") as in_handle_R1, \
     gzip.open(output_file_R1, "wt") as out_handle_R1, \
     gzip.open(output_file_R2, "wt") as out_handle_R2:
    for title, seq, qual in FastqGeneralIterator(in_handle_R1):
        read_count += 1
        # if read_count == 10000: # Debuging purposes
        #     break

        # Find spacers
        sp1   = find_seq(spacer1, seq)
        sp2   = find_seq(spacer2, seq)
        sp3   = find_seq(spacer3, seq)
        sp_P7 = find_seq(spacer_P7, seq)

        if not sp1:
            # sys.stdout.write('\n')
            statistics['Spacer_1 not found'] += 1
            continue
        if not sp2:
            # sys.stdout.write('\n')
            statistics['Spacer_2 not found'] += 1
            continue
        if not sp3:
            # sys.stdout.write('\n')
            statistics['Spacer_3 not found'] += 1
            continue

        bcd1     = seq[:sp1]
        bcd1_q   = qual[:sp1]

        bcd2     = seq[sp1 + len(spacer1):sp1 + len(spacer1) + 8]
        bcd2_q   = qual[sp1 + len(spacer1):sp1 + len(spacer1) +8]

        bcd3     = seq[sp2 + len(spacer2):sp2 + len(spacer2) + 8]
        bcd3_q   = qual[sp2 + len(spacer2):sp2 + len(spacer2) + 8]

        new_read   = seq[sp3 + spacer3_len:sp_P7]    # If sp_P7 is not found then: sp_P7 = None
        new_read_q = qual[sp3 + spacer3_len:sp_P7]   # If sp_P7 is not found then: sp_P7 = None

        barcode    = bcd2 + bcd1
        barcode_q  = bcd2_q + bcd1_q

        if len(bcd1) != 8 or len(bcd2) != 8 or len(bcd3) != 8:
            statistics['Barcode wrong length'] += 1
            continue

        if Levenshtein.distance(bcd3,args.barcode) > 1:
            continue
        statistics['Reads with modality barcode'] += 1

        if len(new_read) < 18:
            statistics['Read too short'] += 1
            continue

        statistics['Good read in modality'] += 1

        out_handle_R1.write("@{}\n{}\n{}\n{}\n".format(title, new_read,'+', new_read_q))
        out_handle_R2.write("@{}\n{}\n{}\n{}\n".format(title, barcode,'+', barcode_q))


    yaml.dump(statistics, sys.stderr)
