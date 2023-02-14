# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 17:33:10 2021

@author: Danrui Wang
@email: wangdanrui18@mails.ucas.ac.cn

"""

import pandas as pd
import sys


def length_calcu(seq_fa):
    refs = []
    l_refs = []
    with open(seq_fa, "r") as f1:
        data = f1.readlines()
        for i in data:
            i = i.replace("\n", "")
            if i.startswith(">"):
                i = i.replace(">","",1)
                refs.append(i)
            else:
                l = len(i)
                l_refs.append(l)
    max_len = max(l_refs)
    l_facts = []
    for i in l_refs:
        l_fac = max_len / i
        l_facts.append(l_fac)
    dict2 = dict(zip(refs, l_facts))
    return dict2



def reads_calcu(in_reads):
    with open(in_reads, "r") as r1:
        samples = []
        n_reads = []
        data = r1.readlines()
        for line in data:
            sample = line.split("\t")[0]
            n_read = line.split("\t")[-1].replace("\n", "")
            samples.append(sample)
            n_reads.append(float(n_read))
    max_size = max(n_reads)
    s_facts = []
    for i in n_reads:
        s_fac = max_size / i
        s_facts.append(s_fac)
    dict1 = dict(zip(samples, s_facts))
    return samples, dict1

def normal_abun(abun_in, read_in, length_in, normal_out):
    dict2 = length_calcu(length_in)
    samples, dict1 = reads_calcu(read_in)
    abun = pd.read_csv(abun_in, sep="\t", index_col=0)
    col_index = abun.index.values
    col_names = abun.columns.values
    normal_df = pd.DataFrame(index=abun.index, columns=abun.columns)
    def normal_reads(sam):
        for i in col_index:
            normal_df.loc[i, sam] = abun.loc[i, sam] * dict1[sam]
    def normal_length(ref):
        for i in col_names:
            normal_df.loc[ref, i] = abun.loc[ref, i] * dict2[ref]
    for sam in samples:
        normal_reads(sam)
    for ref in col_index:
        normal_length(ref)
    normal_df.to_csv(normal_out, sep="\t", index=True)

def int_abun(normal_in, int_out):
    normal_df = pd.read_csv(normal_in, sep="\t", index_col=0)
    int_df = normal_df.round()
    int_df.to_csv(int_out, sep="\t", index=True)



def main(argv=sys.argv[1:]):
    f_abun = sys.argv[1]
    f_read = sys.argv[2]
    f_leng = sys.argv[3]
    f_normal = sys.argv[4]
    f_round = sys.argv[5]
    normal_abun(f_abun, f_read, f_leng, f_normal)
    int_abun(f_normal,f_round)

main()