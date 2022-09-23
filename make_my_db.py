#!/usr/bin/env python3
import sys
import os
import subprocess
import re
import random
from optparse import OptionParser

USAGE = """usage: %prog [OPTIONS]

%prog creates a reference database and the necessay indices for use by
EMIRGE from your reference database. 
Requires bowtie2-build can be found in path

"""
def INFO(information):
    print(information)
    return

def pairs(lst):
    """Creates pairwise iterator on `lst`.
    I.e. if lst returns 1..5, pairs(lst) will return (1,2), (3,4), (5,"")
    """
    it = iter(lst)
    for item in it:
        try:
            yield item, next(it)
        except StopIteration:
            yield item, ""

def cluster_fasta(vsearch_bin, filein, minlen, maxlen, clusterid, threads):
    """Runs vsearch cluster_fast on `filein`, considering only sequences
    with `minlen` <= sequence length <= `maxlen`
    """
    filein_split = filein.split(".")
    name_split = filein_split[0].split("_")
    fileout = name_split[-1] + "." + ".".join(
        ["ge{0}bp".format(minlen), "le{0}bp".format(maxlen), str(clusterid)] +
        filein_split[-2:-1]
    )

    if os.path.isfile(fileout) and os.path.getsize(fileout) >0:
        print ("Found existing file \"{0}\". Skipping clustering stage."
               .format(fileout))
        return fileout

    cmd = [vsearch_bin,
           "--threads", str(threads),
           "--minseqlength", str(minlen),
           "--maxseqlength", str(maxlen),
           "--fasta_width", "0",
           "--notrunclabels",
           "--centroids", fileout,
           "--cluster_fast", filein,
           "--id", str(clusterid)]

    print (" ".join(["Running: "] + cmd))
    subprocess.call(cmd)
    print("")

    return fileout

def randomize_ambiguous(seq):
    """Replaces IUPAC ambiguous characters in seq with random characters,
    respecting options, i.e. "R" replaced with A or G.
    """
    iupac_map = {
    "": " ",
    "R": "AG",    # Purine (A or G)
    "Y": "CT",    # Pyrimidine (C, T, or U)
    "M": "CA",    # C or A
    "K": "TG",    # T, U, or G
    "W": "TA",    # T, U, or A
    "S": "CG",    # C or G
    "B": "CTG",   # C, T, U, or G (not A)
    "D": "ATG",   # A, T, U, or G (not C)
    "H": "ATC",   # A, T, U, or C (not G)
    "V": "ACG",   # A, C, or G (not T, not U)
    "N": "ACTG"   # Any base (A, C, G, T, or U)
    }

    # slower (but clearer) implementation:
    # return "".join([choice(iupac_map.get(x, "ACTG")) for x in seq.upper()])
    seq = seq.upper().replace("U", "T")
    list = [ok + random.choice(iupac_map.get(fix, "ACGT"))
            for ok, fix in pairs(re.split("([^ACGT])", seq))]
    return "".join(list).rstrip(" ")

def randomize_ambiguous_fasta(filein, folder=None):
    """Replaces IUPAC ambiguous character in `filein` with random choice
    from allowed characters.
    Returns output filename.
    """
    filein_split = filein.split(".")
    fileout = ".".join(filein_split[0:-1] + ["fixed"] + filein_split[-1:])
    if folder is not None:
        fileout = os.path.join(folder, os.path.basename(fileout))

    if os.path.exists(fileout) and os.path.getsize(fileout) > 0:
        print ("Found existing file \"{0}\". Skipping randomization stage."
               .format(fileout))
        return fileout

    processed = 0
    total = os.path.getsize(filein)
    dots = 0
    linewidth = 77
    print("Randomizing ambiguous bases")
    print("|" + "-" * (linewidth-2) + "|")
    with open(filein, "r") as inf, open(fileout, "w") as outf:
        for line in inf:
            if line[0] == '>':
                outf.write(line)
            else:
                outf.write(randomize_ambiguous(line.rstrip("\n")))
                outf.write("\n")
            processed += len(line)
            dotstoprint = int(linewidth * processed / total) - dots
            sys.stderr.write("."*dotstoprint)
            sys.stderr.flush()
            dots += dotstoprint
    sys.stderr.write("\n")
    return fileout

#修改成调用bowtie2构建索引
def build_bowtie2_index(bowtie2_bin, filein):
    """Calls bowtie-build on `filein` to compute bowtie index"""
    fileout = ".".join(filein.split(".")[:-1])
    cmd = [bowtie2_bin, "-o", "0", filein, fileout]
    print("Running: " + " ".join(cmd))
    subprocess.call(cmd)


def main(argv=sys.argv[1:]):
    parser = OptionParser(USAGE)
    parser.add_option(
        "-p", "--threads", dest="threads", type="int",
        default=0,
        help="number of threads to use for vsearch clustering of database (default = use all available)"
    )
    parser.add_option(
        "-t", "--tmpdir", dest="tmpdir", metavar="DIR", type="string",
        default="/tmp",
        help="working directory for temporary files (default = %default)"
    )
    parser.add_option(
        "-m", "--min-len", dest="min_len", metavar="LEN", type="int",
        default=1200,
        help="minimum reference sequence length, no less than 1024 (default = %default)"
    )
    parser.add_option(
        "-M", "--max-len", dest="max_len", metavar="LEN", type="int",
        default=2000,
        help="maximum reference sequence length (default = %default)"
    )
    parser.add_option(
        "-i", "--id", dest="clusterid", metavar="FLOAT", type="float",
        default=0.97,
        help="Cluster at this fractional identity level (default = %default)"
    )
    parser.add_option(
        "-k", "--keep", dest="keep", action="store_true",
        help="keep intermediary files (default: do not keep)"
    )
    parser.add_option(
        "-V", "--vsearch", dest="vsearch", metavar="FILE", type="string",
        default="vsearch",
        help="path to vsearch binary (default: look in $PATH)"
    )
    parser.add_option(
        "-B", "--bowtie2-build", dest="bowtie2", metavar="FILE", type="string",
        default="bowtie2-build",
        help="path to bowtie2-build binary (default: look in $PATH)"
    )
    parser.add_option(
        "-d", "--my-database", dest="mydb", type="string",
        help="input your reference database of functional gene"
    )
    (options, args) = parser.parse_args(argv)
    clustered_fasta = cluster_fasta(options.vsearch, options.mydb,
                                    options.min_len, options.max_len,
                                    options.clusterid, options.threads)
    randomized_fasta = randomize_ambiguous_fasta(clustered_fasta,
                                                 folder="")
    build_bowtie2_index(options.bowtie2, randomized_fasta)


if __name__ == '__main__':
    main()

