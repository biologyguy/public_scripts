#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, version 3 of the License (GPLv3).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details at http://www.gnu.org/licenses/.

Copyright (C) Steve Bond, 2014

name: siRNA_predict.py
date: Nov-27-2014
version: 1.0
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/public_scripts
© license: Gnu General Public License, Version 3.0 (http://www.gnu.org/licenses/gpl.html)
derivative work: No

Description:
Python implementation of siRNA prediction criteria as determined by Reynolds et al., 2004, Nat Biotechnol 22(3):326-330.
This script will generate equivalent output to the online tool at http://tiny.cc/naus_siRNA_pred. The sense strand cDNA
is provided as input on the command line, either as a string or FASTA formatted file, and the output is a list of scores
for all possible 19-mer sequences, as a preformatted table or in CSV format. For a detailed description of the
parameters the script takes, navigate to the directory containing the program within a terminal window, and run the
following command:

    python ./siRNA_predict.py -h

"""
from re import findall, sub
import argparse
from os.path import isfile

parser = argparse.ArgumentParser(prog="siRNA prediction",
                                 description="Implementation of siRNA design algorithm developed by "
                                             "Reynolds et al., 2004, Nat Biotechnol 22(3):326-330",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('sequence', help='Input DNA sequence to analyze (FASTA file or string)')
parser.add_argument('-c', '--csv', help='Output as pure CSV', action="store_true")

in_args = parser.parse_args()


def si_score(sequence):
    _score = 0
    # Criteria 1: Moderate to low (30%-52%) GC Content -> 1 point
    gc_count = len(findall("G", sequence)) + len(findall("C", sequence))
    if 6 <= gc_count <= 10:
        _score += 1
    
    # Criteria 2: At least 3 A/Us at positions 15-19 (sense) -> 1 point per (A/U)
    terminal_at_count = len(findall("A", sequence[14:])) + len(findall("T", sequence[14:]))
    _score += terminal_at_count
    
    # Criteria 3: Contains stretch of 4 or more bases, such as AAAA or CCCC -> -1
    if len(findall("AAAA", sequence)) > 0 or len(findall("TTTT", sequence)) > 0 or len(findall("GGGG", sequence)) > 0\
            or len(findall("CCCC", sequence)) > 0:
        _score -= 1
    
    # Criteria 4: A at position 19 (sense) -> 1 point
    if sequence[18] == "A":
        _score += 1
        
    # Criteria 5: A at position 3 (sense) -> 1 point
    if sequence[2] == "A":
        _score += 1
    
    # Criteria 6: T at position 10 (sense) -> 1 point
    if sequence[9] == "T":
        _score += 1
        
    # Criteria 7: G/C at position 19 (sense) -> -1 point
    if len(findall("[GC]", sequence[18])) > 0:
        _score -= 1

    # Criteria 8: G at position 13 (sense) -> -1 point
    if sequence[12] == "G":
        _score -= 1
    return _score


if isfile(in_args.sequence):
    with open(in_args.sequence, "r") as ifile:
        full_seq = ifile.read()

    full_seq = sub(">.*", "", full_seq).upper()
    full_seq = sub("\s[0-9]", "", full_seq)
    full_seq = sub("U", "T", full_seq)
    full_seq = sub("[^ATCG]", "X", full_seq)

else:
    full_seq = in_args.sequence.upper()
    full_seq = sub("U", "T", full_seq)
    full_seq = sub("[^ATCG]", "X", full_seq)

si_seqs_list = [[], [], [], [], [], [], [], [], [], []]
for i in range(len(full_seq) - 18):
    seq = full_seq[i:i + 19]
    score = si_score(seq)
    score = score if score >= 0 else 0
    si_seqs_list[score].append((seq, i + 1))

biggest_column = 0
for si_seq in si_seqs_list:
    if len(si_seq) > biggest_column:
        biggest_column = len(si_seq)

# super clunky text formating... But looks good in terminal
output = "Score:\t\t9\t\t\t\t8\t\t\t\t7\t\t\t\t6\t\t\t\t5\t\t\t\t4\t\t\t\t3\t\t\t\t2\t\t\t\t1\t\t\t\t0\n"
for row in range(biggest_column):
    for index in [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]:
        if row < len(si_seqs_list[index]):
            output += "%s— %s\t" % (str(si_seqs_list[index][row][1]).ljust(5), si_seqs_list[index][row][0])
        else:
            output += "\t\t\t\t"
    output += "\n"

if in_args.csv:
    output = sub("\t{4}", ",,", output)
    output = sub("\t{2}", ",", output)
    output = sub("\t", ",", output)
    output = sub("—", ",", output)
    output = sub(" ", "", output)

print(output)
