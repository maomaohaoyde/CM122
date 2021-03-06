import pandas as pd
import os
import sys
from collections import defaultdict, Counter
from os.path import join,exists,splitext
import math
import numpy as np
import itertools
def read_file(read_fn):
    # helper funtion to read both reads and strains files
    f = open(read_fn, 'r')
    first_line = True
    all_reads = []
    for line in f:
        if first_line:
            first_line = False
            continue
        line = line.strip()
        all_reads.append(line)
    return all_reads


def find_SNPs(input_strains):
    # compare four strains and find the positions of snps
    snps = []
    size = len(input_strains[0])
    index = 0
    for i in range(size):
        if input_strains[0][i] == input_strains[1][i] and input_strains[2][i] == input_strains[3][i] and input_strains[0][i] == input_strains[2][i]:
            continue
        else:
            snps.append(i)
    return snps




def alternative_hash(snps, input_strains):
    # hash substrings of strains that could contain the snps
    strain_hash = defaultdict(tuple)
    for j in range(len(input_strains)):
        for key in snps:
            for i in range(key-49,key+1):
                ref_piece = input_strains[j][i:i+50]
                strain_hash[ref_piece]=(i,key-i)
    return strain_hash

#   recursively calculate the max determinant, is not used in the method
# def max_det(strain_matrix,left,flag,snps_freq):
#     if left == 0:
#         hello = [tuple(strain) for strain in strain_matrix]
#         hi = np.array(hello)
#         bye = np.dot(hi.T,hi)
#         return np.linalg.det(bye) , strain_matrix, snps_freq
#     else:
#         if flag:
#             k = left
#             snps_freq[k] = 1 - snps_freq[k]
#             for i in range(len(strain_matrix[k])):
#                 if strain_matrix[k][i] == 1:
#                     strain_matrix[k][i] = 0
#                 else:
#                     strain_matrix[k][i] = 1
#         a , a1, a2 = max_det(strain_matrix,left-1,1,snps_freq)
#         b , b1, b2 = max_det(strain_matrix,left-1,0,snps_freq)
#         if a > b:
#             c = a
#             d = a1
#             e = a2
#         else:
#             c = b
#             d = b1
#             e = b2
#         return c,d,e

def sum_freq(strain_matrix,left,flag,snps_freq):
    # recursively calculate the sum of the frequency calculated from the different combination of distribution of snps
    if left == 0:
        strains_snps = [tuple(d1) for d1 in strain_matrix]
        strains = np.array(strains_snps)
        snp_freqs1 = np.array(snps_freq)


        strain_freqs = np.linalg.lstsq(strains, snp_freqs1)

        for i in range(len(strain_freqs[0])):
            if strain_freqs[0][i] < 0:
                strain_freqs[0][i] = 0 - strain_freqs[0][i]
        return strain_freqs[0]
    else:
        if flag:
            k = left
            snps_freq[k] = 1 - snps_freq[k]
            for i in range(len(strain_matrix[k])):
                if strain_matrix[k][i] == 1:
                    strain_matrix[k][i] = 0
                else:
                    strain_matrix[k][i] = 1
        a = sum_freq(strain_matrix,left-1,1,snps_freq)
        b = sum_freq(strain_matrix,left-1,0,snps_freq)
        c = []
        for i in range(len(a)):
            c.append(a[i]+b[i])

        return c

def read_mapping(ht, reads, input_strains,snps):
    mapping  = defaultdict(list)
    strain_snps = [[] for i in range(len(snps))]
    input_counts = defaultdict(list)

    #find the characters at the position of snps
    for i in range(len(reads)):
        if reads[i] in ht:
            a = ht[reads[i]][0]
            b = ht[reads[i]][1]
            key = a+b
            mapping[key].append(reads[i][b])

    counter = 0

    #assume for the moment the first strain has all the snps
    #update if other strains contain snps accordingly
    #record every appearance of the character at the position of snps in the first strain
    for k in snps:
        cnt = Counter()
        for word in mapping[k]:
            cnt[word] += 1
        temp = input_strains[0][k]
        strain_snps[counter].append(1)
        for i in range(1,len(input_strains)):
            if input_strains[i][k] == temp:
                strain_snps[counter].append(1)
            else:
                strain_snps[counter].append(0)

        for j in range(len(mapping[k])):
            if mapping[k][j] == temp:
                input_counts[counter].append(1)
            else:
                input_counts[counter].append(0)
        counter += 1


    #calculate the snps frequency array
    snp_freqs = []
    for k in range(len(snps)):
        raw_snps = input_counts[k]
        snp_freq = float(sum(raw_snps))/len(raw_snps)
        snp_freqs.append(snp_freq)


    # calculate the sum of calculated strain frequency with sum_freq
    # here the average of them is used as the final result
    a = sum_freq(strain_snps, 6, 1, snp_freqs)
    b = sum_freq(strain_snps, 6, 0, snp_freqs)
    c = []
    for i in range(len(a)):
        c.append((a[i] + b[i])/128)

    return c



if __name__ == "__main__":
    input_folder = './hw4_W_1/'
    f_base = 'hw4_W_1'
    read_fn_end = '{}_reads.txt'.format(f_base)
    read_fn = join(input_folder,read_fn_end)
    strain_fn_end = '{}_strains.txt'.format(f_base)
    strain_fn= join(input_folder,strain_fn_end)
    strains = read_file(strain_fn)
    reads = read_file(read_fn)
    snps = find_SNPs(strains)
    hashing = alternative_hash(snps,strains)
    strain_freqs = read_mapping(hashing,reads,strains,snps)
    output_result = ''
    for i in range(len(strains)):
        output_result = output_result + str(strain_freqs[i]) + ',' + strains[i] + '\n'
    output_fn = '{}_ans.txt'.format(f_base)
    output_str = '>hw4_W_1'+'\n'+output_result
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)

