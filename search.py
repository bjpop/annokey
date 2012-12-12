#!/bin/env python

import sys
import csv
from argparse import ArgumentParser
import re
'''
Version modified from previous to accept Entrez_parser.py parsed xml file such as resulting from Homo_sapiens.xml. Accepts this output file, searchKeys and genes in tsv format. Works offline so that any number of variants can be tested.
'''

def parse_args():
    parser = ArgumentParser(
        description = 'Search NCBI for genes of interest, based on concept-keyword search')
    parser.add_argument(
        '--searchKeys', metavar = 'SEARCHKEYS', type = file, required=True,
        help ='The tab separated file containing name of the keywords to be searched')
    parser.add_argument(
        '--genes', metavar = 'GENES', type = file, required=True,
        help ='The tab separated file containing the gene information including name of the gene, one gene name per line')
    parser.add_argument(
        '--geneDb', metavar = 'GENEDB', type = file, required=True,
        help ='The tab separated file containing the gene information provided by ncbi, parsed representation of Hsxml')
    return parser.parse_args()


def read_searchKeys(file):
    searchItems = []
    for line in file:
        term = line.strip()
        searchItems.append(term)
    return searchItems


def search_genes(genefile, keywords, Hsdict):
    '''Searches gene information from HsDict, and
       prints the results
       1. searches keywords, and appends results to the original data
          ['Keyword', 'Key Rank']
    '''
    additionalFields = ['Keyword', 'Key Rank']

    genes = csv.DictReader(genefile, delimiter='\t')
    
    # Gets headers and 
    # inserts addtitional headers for search results
    # at the end of headers
    headers = genes.fieldnames
    newHeaders = headers + additionalFields
    
    # Prints headings 
    # (Finally, will write a result file instead of printing)
    print '\t'.join(newHeaders)

    # Searches gene information from gene database
    for gene in genes:
        geneInfo = [gene[field] for field in headers]
     
        if gene['Func'] == 'intergenic':
            # Skip
            geneInfo += ['', '']
            continue
        else:
            # Search
            try:
                geneName = gene['Gene']
                dbEntry = Hsdict[geneName]
                searchResults = search_keywords(dbEntry, keywords)
                if searchResults != None:
                    keyword, rank = searchResults
                    geneInfo += [keyword, str(rank)]
                else:
                    # if db does not contain the keywords
                    geneInfo += ['Not Found', '']
            except KeyError:
                # if db does not contain the geneName
                geneInfo += ['Gene Not Found', '']
     
        print '\t'.join(geneInfo)


    

def search_keywords(dbEntry, keywords):
    '''Searches keyword in the order (according to the rank) and
       returns the keyword and its rank (position in the list)
       If fail to search, (python) returns None.
    '''
    for n,item in enumerate(keywords):
        if item in dbEntry:
            return (item, n+1)
    


def read_geneDb_asDictionary(geneDb):
    '''Creates gene dictionary from gene database'''
    geneDict = {}
    for line in geneDb:
        linelist = line.split('\t')
        # .tsv file has no heading
        # index 1: gene name
        # there is a case linelist is empty..(so, try and except) 
        try:
            genename = linelist[1]
            geneDict[genename] = line
        except:
            continue

    return geneDict


def main():
    args = parse_args()
    geneDict = read_geneDb_asDictionary(args.geneDb)
    keywords = read_searchKeys(args.searchKeys)
    search_genes(args.genes, keywords, geneDict)

if __name__ == '__main__':
    main()


