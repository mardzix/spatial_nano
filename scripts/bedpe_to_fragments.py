#!/usr/bin/env python3

# Use bedpe as input, generate with such command:
# bedtools bamtobed -bedpe -i {file} 2>/dev/null | python3 bedpe_to_fragments.py

import sys

for line in sys.stdin:
    line = line.rstrip().split('\t')
    start1 = int(line[1])
    start2 = int(line[4])
    end1 = int(line[2])
    end2 = int(line[5])
    chr1 = line[0]
    chr2 = line[3]
    name = line[6]
    score = line[7]
    strand = '.'
    if chr1 != chr2:
        continue
    new_start = min(start1, start2, end1, end2)
    new_end   = max(start1, start2, end1, end2)
    sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr1, new_start,new_end,name,score,strand))
