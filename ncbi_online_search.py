#!/bin/env python

import sys
import csv
import itertools as iter
from argparse import ArgumentParser
from Bio import Entrez
import logging
import pickle
import pprint

Entrez.email = "bjpope@unimelb.edu.au"

def readSearchKeys(filename):
    searchItems = {}
    with open(filename) as file:
        reader = csv.reader(file, delimiter=',', quotechar='"')
        header = reader.next()
        for h in header:
            searchItems[h] = []
        for line in reader:
            for n,item in enumerate(line):
                if item:
                    searchItems[header[n]].append(item)
    return searchItems

def readGenes(filename):
    result = []
    with open(filename) as file:
        for line in file:
            result.append(line.strip())
    return result

def getGeneIDs(genes):
    resultIDs = []
    for g in genes:
        searchTerm = "Homo sapiens[organism] AND {}[sym]".format(g)
        request = Entrez.esearch("gene", term=searchTerm)
        result = Entrez.read(request)
        resultIDs.extend(result['IdList'])
    return resultIDs

def searchGenes(ids):
    idString = ','.join(ids)
    request = Entrez.epost("gene",id=idString)
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.efetch(db="gene", webenv=webEnv, query_key=queryKey, retmode='xml')
    #data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey, retmode='xml')
    allResults = Entrez.read(data)
    # return filter(isHumanResult, allResults)
    return allResults

def isHumanResult(result):
    return result['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_taxname'] == 'Homo sapiens'

def parseArgs():
    parser = ArgumentParser(
        description = 'Search NCBI for genes of interest, based on concept-keyword search')
    parser.add_argument(
        '--searchKeys', metavar = 'SEARCHKEYS', type = str,
        required=True, help ='name of the search keywords CSV file')
    parser.add_argument(
        '--genes', metavar = 'GENES', type = str,
        required=True, help ='name of the gene list text file, one gene name per line')
    parser.add_argument(
        '--saveCache', metavar = 'CACHEFILE', type = str,
        help ='Save a cache of the downloaded results from NCBI into this file')
    parser.add_argument(
        '--loadCache', metavar = 'CACHEFILE', type = str,
        help ='Load previously cached download results from NCBI from this file')
    parser.add_argument('--logFile', metavar='FILENAME', type=str,
                    help='log progress in FILENAME')
    return parser.parse_args()

def main():
    args = parseArgs()
    if args.logFile:
        logging.basicConfig(filename=args.logFile,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {}'.format(' '.join(sys.argv)))
    searchKeys = readSearchKeys(args.searchKeys)
    logging.info('{} search keys read from {}\n{}'.format(len(searchKeys),
       args.searchKeys, ','.join(searchKeys.keys())))
    geneNames = readGenes(args.genes)
    logging.info('{} genes read from {}'.format(len(geneNames), args.genes))
    if args.loadCache:
        with open(args.loadCache) as file:
            searchResults = pickle.load(file)
        logging.info('{} search results loaded from cache file {}'.format(len(searchResults), args.loadCache))
    else:
        geneIds = getGeneIDs(geneNames)
        logging.info('{} gene IDs retrieved from NCBI\n{}'.format(len(geneIds),
            ','.join(geneIds)))
        searchResults = searchGenes(geneIds) 
        logging.info('{} search results retrieved from NCBI'.format(len(searchResults)))
        if args.saveCache:
            with open(args.saveCache, 'w') as file:
                pickle.dump(searchResults, file)
            logging.info('{} search results saved to cache file {}'.format(len(searchResults), args.saveCache))
    for result in searchResults:
        print(result)

if __name__ == '__main__':
    main()
