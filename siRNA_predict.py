# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details at http://www.gnu.org/licenses/.
#
# Copyright (C) Steve Bond, 2014
# Contact: biologyguy@gmail.com

"""
Description:
Python implementation of siRNA prediction criteria as determined by Reynolds et al., 2004, Nat Biotechnol 22(3):326-330.
This script will generate equivalent output to the online tool at http://tiny.cc/naus_siRNA_pred. The sense strand cDNA
is provided as input on the command line, either as a string or FASTA formatted file, and the output is scores for all
possible 19-mer sequences, as a preformatted table or in CSV format. For a detailed description of the parameters the
script takes, navigate to the directory containing the program within a terminal window, and run the following command
(tested using python 2.7.8 on Windows 7, Linux Mint 17, and Mac OS X 10.9):

python ./siRNA_predict.py -h

"""
import re
import sys
import argparse

parser = argparse.ArgumentParser(prog="siRNA prediction",
                                 description="Implementation of siRNA design algorithm developed by Reynolds et al., 2004, Nat Biotechnol 22(3):326-330",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-s', '--sequence', metavar='', help='Input DNA sequence to analyze', default=False)
parser.add_argument('-c', '--csv', help='Output as pure CSV', action="store_true", default=False)
parser.add_argument('-f', '--fasta', metavar='', help='Read in sequence, either raw or FASTA', action="store", default=False)

incoming_args = parser.parse_args()

if incoming_args.sequence:
    full_seq = incoming_args.sequence.upper()
    full_seq = re.sub("[^ATCG]", "X", full_seq)

elif incoming_args.fasta:
    with open(incoming_args.fasta, "r") as ifile:
        full_seq = ifile.read()

    full_seq = re.sub(">.*\n", "", full_seq).upper()
    full_seq = re.sub("[^ATCG]", "X", full_seq)

else:
    sys.exit("You need to provide a sequence, using either the -s or -f flag.")


def si_score(sequence):
    score = 0
    #Criteria 1: Moderate to low (30%-52%) GC Content -> 1 point
    GC_count = len(re.findall("G", sequence)) + len(re.findall("C", sequence))
    if 6 <= GC_count <= 10:
        score += 1
    
    #Criteria 2: At least 3 A/Us at positions 15-19 (sense) -> 1 point per (A/U)
    terminal_AT_count = len(re.findall("A", sequence[14:])) + len(re.findall("T", sequence[14:]))
    score += terminal_AT_count
    
    #Criteria 3: Contains stretch of 4 or more bases, such as AAAA or CCCC -> -1
    if len(re.findall("AAAA", sequence)) > 0 or len(re.findall("TTTT", sequence)) > 0 or len(re.findall("GGGG", sequence)) > 0 or len(re.findall("CCCC", sequence)) > 0:
        score -= 1       
    
    #Criteria 4: A at position 19 (sense) -> 1 point
    if sequence[18] == "A":
        score += 1
        
    #Criteria 5: A at position 3 (sense) -> 1 point
    if sequence[2] == "A":
        score += 1    
    
    #Criteria 6: T at position 10 (sense) -> 1 point
    if sequence[9] == "T":
        score += 1
        
    #Criteria 7: G/C at position 19 (sense) -> -1 point
    if len(re.findall("[GC]", sequence[18])) > 0:
        score -= 1

    #Criteria 8: G at position 13 (sense) -> -1 point
    if sequence[12] == "G":
        score -= 1

    return score


si_seqs_list = [[], [], [], [], [], [], [], [], [], []]
for next in range(len(full_seq) - 18):
    seq = full_seq[next:next + 19]
    score = si_score(seq)
    score = score if score >= 0 else 0
    si_seqs_list[score].append((seq, next + 1))

biggest_column = 0
for next in si_seqs_list:
    if len(next) > biggest_column:
        biggest_column = len(next)

output = "\t\t9\t\t\t\t8\t\t\t\t7\t\t\t\t6\t\t\t\t5\t\t\t\t4\t\t\t\t3\t\t\t\t2\t\t\t\t1\t\t\t\t0\n"
for row in range(biggest_column):
    for index in [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]:
        if row < len(si_seqs_list[index]):
            output += "%s-%s\t" % (str(si_seqs_list[index][row][1]).ljust(5), si_seqs_list[index][row][0])
        else:
            output += "\t\t\t\t"
    output += "\n"

if incoming_args.csv:
    output = re.sub("\t{4}", ",,", output)
    output = re.sub("\t{2}", ",", output)
    output = re.sub("\t", ",", output)
    output = re.sub("-", ",", output)

print(output)
