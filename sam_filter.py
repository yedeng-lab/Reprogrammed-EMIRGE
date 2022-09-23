#!/usr/bin/env python3
import sys
# from optparse import OptionParser
import re

# USAGE = """usage: %prog [OPTIONS]
#
# %prog creates a filtered sam file for use by
# EMIRGE.
#
# """


with open(sys.argv[1], mode='r', encoding='utf-8') as input:
    with open(sys.argv[2], mode='w', encoding='utf-8') as output:
        sys.stdout = output
        for line in input:
            if line.startswith('@'):
                output.write(line)
            else:
                line_cigar = line.split('\t')[5]
                if ('I' in line_cigar) or ('D' in line_cigar) or ('S' in line_cigar) or ('H' in line_cigar):
                    pass

                else:
                    output.write(line)
