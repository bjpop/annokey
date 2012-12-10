#!/bin/env python

import sys
import csv
from argparse import ArgumentParser
from Bio import Entrez
import re

Entrez.email = "djpark@unimelb.edu.au"


def parseArgs():
    parser = ArgumentParser(
        description = 'Search NCBI for genes of interest, based on concept-keyword search')
    parser.add_argument(
        '--searchKeys', metavar = 'SEARCHKEYS', type = file, required=True,
        help ='name of the search keywords CSV file')
    parser.add_argument(
        '--genes', metavar = 'GENES', required=True,
        help ='name of the gene list text file, one gene name per line')
    return parser.parse_args()


def readSearchKeys(file):
    searchItems = []
    file.next()
    for line in file:
        term = str(line.strip())
        searchItems.append(term)
    return searchItems


def readGenes(genes, searchKeys):
    reader = csv.reader(open(genes, 'rU'), delimiter='\t', quotechar='"')
    header = reader.next()
    for h in header:
        if h == 'Gene':
            gene = header.index(h)
        if h == 'Chr':
            chr = header.index(h)
# output the header with additional columns
    header.insert(2,'Keyword')
    header.insert(3,'Key Rank')
    print "\t".join(str(thing) for thing in header)


    for line in reader:
        uid = []
        g = line[gene]
        c = line[chr]
        c = re.sub('chr','',c)
# append line list with two slots: default no keyword hit and ranking
# skip 'intergenic'
        if line[0] == 'intergenic':
            continue
# or carry on happily
        else:
            line.insert(2,'no keyword hit')
            line.insert(3,9999)
            try:
                request = Entrez.esearch(db="gene", term="Homo sapiens[ORGN] AND %s[SYM] AND %s[CHR]" % (g,c), retmax=1)
                result = Entrez.read(request)
                uid.extend(result['IdList'])
                genedbase = searchGene(uid)
                gbasestring = str(genedbase)
                key_in_dbase(gbasestring,searchKeys,g,c,line)
            except:
# append with message that nothing found and return 9999
                line[2] = 'no gene found'
                line[3] = 9999
                pass
            print "\t".join(str(thing) for thing in line)


def searchGene(uid):
    request = Entrez.epost("gene",id=uid)
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.efetch(db="gene", webenv=webEnv, query_key=queryKey, retmode='xml')
    Result = Entrez.read(data)
    return Result


def key_in_dbase(gbasestring,searchKeys,g,c,line):
    for item in searchKeys:
        if item in gbasestring:
                v = searchKeys.index(item)
# append line and output
                line[2] = item
                line[3] = v
                return line
                break    
        else:
            continue
        break


def main():
    args = parseArgs()
    searchKeys = readSearchKeys(args.searchKeys)
    parseGenes = readGenes(args.genes, searchKeys)

if __name__ == '__main__':
    main()


