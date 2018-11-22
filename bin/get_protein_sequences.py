#!/usr/bin/env python

"""
This script downloads desired protein sequences from a given list of taxon IDs.

With the aid of the UniProt API, the sequences of interest are downloaded. The protein of interest is passed by the
parameter `--gene_name`, the organisms of interest are selected based on a text file that contains the tax IDs
(parameter `--input_file`). Optionally, the parameter `--reviewed` can be passed to only download Swiss-Prot entries.
The target for the downloaded file can be passed via `--output_file`.
"""

import argparse
import urllib.parse
import urllib.request

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", action="store", dest="tax_list_file", required=True,
                        help="input file (tax ids, # lines are ignored)")
    parser.add_argument("-o", "--output_file", action="store", dest="fasta_file", required=True,
                        help="output file (fasta file with protein sequences of interest)")
    parser.add_argument("-g", "--gene_name", action="store", dest="gene_name", required=True,
                        help="gene name of interest")
    parser.add_argument("-r", "--reviewed", action="store_true", dest="reviewed", default=False,
                        help="get reviewed entries (Swiss-Prot instead of TrEMBL)")
    args = parser.parse_args()

    output = args.fasta_file

    gene_name = args.gene_name
    reviewed = args.reviewed
    tax_list = []
    tax_list_file = args.tax_list_file
    with open(tax_list_file) as tax_list_open:
        for line in tax_list_open:
            if not line.startswith("#"):
                tax_list.append(line.rstrip())

    taxon_queries = ["taxonomy:%s" % tid for tid in tax_list]
    taxon_query = "(" + " OR ".join(taxon_queries) + ")"
    rev = " AND reviewed:%s" % reviewed if reviewed else ""
    gene_name = " gene_exact:%s" % gene_name
    gene_name_query = " AND " + gene_name
    

    url = "https://www.uniprot.org/uniprot/"
    query = "%s%s%s" % (taxon_query, rev, gene_name_query)
    # print(query)
    params = {"query": query, "force": "yes", "format": "fasta"}
    data = urllib.parse.urlencode(params).encode("utf-8")
    msg = urllib.request.urlretrieve(url=url, filename=output, data=data)[1]
    headers = {j[0]: j[1].strip() for j in [i.split(":", 1)
                                                for i in str(msg).strip().splitlines()]}

if __name__ == "__main__":
    # execute only if run as a script
    main()

