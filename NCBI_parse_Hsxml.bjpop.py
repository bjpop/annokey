#!/bin/env python

import sys
from argparse import ArgumentParser
from Bio import Entrez

Entrez.email = 'djp@unimelb.edu.au'

def parseArgs():
    parser = ArgumentParser(
        description = 'Convert Entrez Gene Homo_sapiens.xml to python dictionary representation')
    parser.add_argument(
        '--Hsxml', metavar = 'HSXML', type = file, required = True,
        help = 'Name of Homo_sapiens.xml file - include a date reference for download for example')
    return parser.parse_args()

args = parseArgs()

#records = Entrez.read(args.Hsxml)
records = Entrez.parse(args.Hsxml)

for record in records:
    genename = record['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
    print genename




        




#def parsegenerecord(gene,chr,EparseHsxml):
##    records = Entrez.parse(Hsxml)
#    for record in EparseHsxml:
#        content = []
#        try:
#            if record[][][][][] == gene:
#                if record[][][]][][] == chr:
#                    #parse contents of record
#                    sym = record[][][][][][][]
#                    chr = record[][][][]][]][]
#                    name = record[][][][][][]
#                    summary = record[][][][][][]
#                    function = record[][][][][][][][]
#                    intns = record[][][][][][][]
#                    generif = record[][][][][][][]
#                    paths = record[][][][][][]
#                    # etc.
#                    # need to concatenate list contents for any of these?
#                    content = [sym,chr,name,summary,function,intns,generif,paths] # etc.
#                    content = "".join(str(thing) for thing in content)
#                    return content
#                else:
#                    continue            
#
#            else:
#                continue
#
#        except:
#            return None
        


#    contents = []
#    gene_name = # record[][][]][][]
#    chr = # record[][][][][][]
#    genchr = str(gene_name)+str(chr)
#    name = # record[][][]][][[]
#    summary = # record[][][[[]][[]][]
#    function = # record[][][][][][][]
#    intns = # record[][][[][][]
#    generif = # record[][]]][][[][][]
#    paths = # record[][][][][][]
#    # and other fields?

 
