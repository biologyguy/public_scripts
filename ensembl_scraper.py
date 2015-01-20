#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, version 2 of the License (GPLv2).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details at http://www.gnu.org/licenses/.

name: siRNA_predict.py
date: Jan-18-2015
version: 1.1
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/public_scripts
© license: Gnu General Public License, Version 2.0 (http://www.gnu.org/licenses/gpl.html)
derivative work: No

Description:
Search ENSEMBL Metazoa for genes and print all returned IDs to a file, sorted by organism. For a detailed description of
the parameters the script takes, navigate to the directory containing the program within a terminal window, and run the
following command:

    python ./ensembl_scraper.py -h
"""

import requests
from bs4 import BeautifulSoup
import argparse
import os
from sys import stdout, exit

parser = argparse.ArgumentParser(prog="ensembl_scraper",
                                 description="Search EnsemblMetazoa for a all genes returned from a search")
parser.add_argument('search_term', help='What would you like to search for?', action='store')
parser.add_argument('-o', '--outfile', help='Send the results to a file, instead of StdOut',
                    action="store", default="%s/ensemble_ids.txt" % os.getcwd())

in_args = parser.parse_args()

# Run the search, and figure out how many pages of results are returned
url = "http://metazoa.ensembl.org/Multi/Search/Results?q=%s;species=all;collection=all;site=ensemblunit"\
      % in_args.search_term
content = requests.get(url).text

soup = BeautifulSoup(content)
try:
    paginate = soup.find('div', {"class": 'paginate'}).find_all('a')
    max_page = 1
    for page in paginate:
        try:
            if int(page.text) > max_page:
                max_page = int(page.text)
        except ValueError:
            continue
except AttributeError:
    max_page = 1

print("%s pages of results were returned" % max_page)
ids = {}
for page_num in range(max_page):
    stdout.write("\rCollecting data from results page %s" % (page_num + 1),)
    stdout.flush()

    url = "http://metazoa.ensembl.org/Multi/Search/Results?page=%s;q=%s;species=all;collection=all;site=ensemblunit"\
          % (page_num + 1, in_args.search_term)
    content = requests.get(url).text

    soup = BeautifulSoup(content)

    for row in soup.find_all('div', {"class": 'row'}):
        sub_soup = BeautifulSoup(str(row))
        lhs = sub_soup.find('div', {"class": "lhs"}).text
        rhs = sub_soup.find('div', {"class": "rhs"}).text

        if lhs == "Gene ID":
            gene_id = rhs

        if lhs == "Species":
            if rhs in ids:
                ids[rhs].append(gene_id)
            else:
                ids[rhs] = [gene_id]

if len(ids) == 0:
    exit("\rNo records found for query '%s'" % in_args.search_term)

output = ""
for species in ids:
    output += "%s\n" % species
    for next_id in ids[species]:
        output += "%s\n" % next_id
    output += "\n"

if in_args.outfile:
    outfile = os.path.abspath(in_args.outfile)
    with open(outfile, "w") as ofile:
        ofile.write(output)
    print("Output written to %s" % outfile)

else:
    print(output)