#!/bin/env python

import sys
import csv
from argparse import ArgumentParser
import re
'''
Version modified from previous to accept Entrez_parser.py parsed xml file such as resulting from Homo_sapiens.xml. Accepts this output file, searchKeys and genes in tsv format. Works offline so that any number of variants can be tested.
'''

def parseArgs():
    parser = ArgumentParser(
        description = 'Search NCBI for genes of interest, based on concept-keyword search')
    parser.add_argument(
        '--searchKeys', metavar = 'SEARCHKEYS', type = file, required=True,
        help ='name of the search keywords CSV file')
    parser.add_argument(
        '--genes', metavar = 'GENES', type = file, required=True,
        help ='name of the gene list text file, one gene name per line')
    parser.add_argument(
        '--Hstsv', metavar = 'HSTSV', type = file, required=True,
        help ='tsv parsed rep of Hsxml')
    return parser.parse_args()


def readSearchKeys(file):
    searchItems = []
    for line in file:
        term = str(line.strip())
        searchItems.append(term)
    return searchItems


def readGenes(genes, searchKeys, Hsdict):
    header = genes.next()
#    header.strip()
    headers = header.split('\t')
    for h in headers:
        h.strip()
        if h == 'Gene':
            gene = headers.index(h)
        if h == 'Chr':
            chr = headers.index(h)
# output the header with additional columns
    headers.insert(2,'Keyword')
    headers.insert(3,'Key Rank')
    print "\t".join(str(thing) for thing in headers)


    for line in genes:
        line = line.split('\t')
        uid = []
        g = line[gene]
        line.insert(2,'no hit')
        line.insert(3,9999)
# append line list with two slots: default no keyword hit and ranking
# skip 'intergenic'
        if line[0] == 'intergenic':
            continue
# or carry on happily
        else:
            try:
                geneentry = str(Hsdict[g])
                key_in_dbase(geneentry,searchKeys,line)
            except:
# append with message that nothing found and return 9999
                line[2] = 'exception'
                line[3] = 9999
            print "\t".join(str(thing) for thing in line)


def key_in_dbase(geneentry,searchKeys,line):
    for item in searchKeys:
        if item in geneentry:
            v = searchKeys.index(item)
# append line and output
            line[2] = item
            line[3] = v
            return line
            break    
        else:
            continue

def HsDict(Hstsv):
    Hsdict = {}
    for line in Hstsv:
        linelist = line.split('\t')
        try:
            genename = linelist[1]
            Hsdict[genename] = line
        except:
            continue
#    print Hsdict
    return Hsdict


def main():
    args = parseArgs()
    Hsdict = HsDict(args.Hstsv)
    searchKeys = readSearchKeys(args.searchKeys)
    parseGenes = readGenes(args.genes, searchKeys, Hsdict)

if __name__ == '__main__':
    main()


