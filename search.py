#!/bin/env python

import sys
import csv
from argparse import ArgumentParser
import re
import pickle
from Bio import Entrez
from lxml import etree

'''
Version modified from previous to accept Entrez_parser.py parsed xml file such as resulting from Homo_sapiens.xml. Accepts this output file, searchKeys and genes in tsv format. Works offline so that any number of variants can be tested.
'''


# for Entrez
Entrez.email = "sorik@student.unimelb.edu.au"

def get_geneIds(genes):
    '''Gets gene Ids through search query.
       The genes are limited to Homo sapiens. 
    '''
    geneIds = []
    for gene in genes:
        searchTerm = 'Homo sapiens[organism] AND {}[sym]'.format(gene['Gene'])
        request = Entrez.esearch(db='gene', term=searchTerm)
        result = Entrez.read(request)
        geneIds.extend(result['IdList'])
    
   
    return geneIds


def fetch_geneInformation(genes):
    '''Searches gene information over online (NCBI).
           genes: csv.DictReader, keywords: list
        1. get gene Ids 
        2. ePost using gene Ids -> get WebEnv and query_key
        3. eFetch using WebEnv and query_key -> get xml containing information
    '''
    
    geneIds = get_geneIds(genes)

    idString = ','.join(geneIds)
    postRequest = Entrez.epost(db='gene', id=idString)
    postResult = Entrez.read(postRequest)
    
    webEnv = postResult['WebEnv']
    queryKey = postResult['QueryKey']
    fetchRequest = Entrez.efetch(db='gene',
                                 webenv=webEnv,
                                 query_key=queryKey,
                                 retmode='xml')
    fetchResults = Entrez.read(fetchRequest)

    return fetchResults



# gene Id = 'Entrezgene_track-info/Gene-track/Gene-track_geneid'
# gene name = 'Entrezgene_gene/Gene-ref/Gene-ref_locus'
# gene synonyms = 'Entrezgene_gene/Gene-ref/Gene-ref_syn'
# gene name(description) = 'Entrezgene_gene/Gene-ref/Gene-ref_desc'
# gene alter name = 'Entrezgene_prot/Prot-ref/Prot-ref_name'
# summary = 'Entrezgene_summary'
# gene RIFs = 'Entrezgene_comments', key==Gene-commentary_text
# function
# component
# process
# pathway = 'Entrezgene_comments/Gene-commentary_comment/Gene-commentary_text'
# interaction
# conserved
# PubMed Id

def get_geneName(geneEntry):
    name = geneEntry.xpath('./Entrezgene_gene/Gene-ref/Gene-ref_locus')
    if len(name) > 0:
        return name[0].text   
  
 

def get_geneContent(geneEntry):
    '''iterchildren'''
    geneContent = {}

    for elem in geneEntry.iterchildren():
        if elem.tag == 'Entrezgene_track-info':
            print elem.tag
        elif elem.tag == 'Entrezgene_gene':
            ref = elem.find('.//Gene-ref')
            name = ref.find('.//Gene-ref_locus')
            desc = ref.find('.//Gene-ref_desc')
            syn = ref.find('.//Gene-ref_syn')
            
        elif elem.tag == 'Entrezgene_summary':
            print elem.tag
        elif elem.tag == 'Entrezgene_comments':
            print elem.tag

        elif elem.tag == 'Entrezgene_prot':
            print elem.tag






def search_genes_xml(genes, keywords, xmlfile):
    '''  '''
     
    genes_new = {}
    content = etree.iterparse(xmlfile, events=('start'), tag='Entrezgene')
    for event, geneEntry in content:
        # get gene name 
        geneName = get_geneName(geneEntry)
        if genes.has_key(geneName):
             gene = genes[geneName]
             if gene['Func'] == 'intergenic':
                 continue
             else: 
                 dbEntry = get_geneContent(geneEntry) 
                 searchResult = search_keywords(geneContent, keywords)
                
        # clear unneeded node references
        geneEntry.clear()
        while geneEntry.getprevious() is not None:
           del geneEntry.getparent()[0]
    del context
        
      
def search_genes(genes, keywords, geneDict):
    '''Searches gene information from geneDict, and
       prints the results.
           genes: csv.DictReader, keywords: list, geneDict: dict
       1. searches keywords, and appends results to the original data
          ['Keyword', 'Key Rank']
    '''
    additionalFields = ['Keyword', 'Key Rank']
    
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
                dbEntry = geneDict[geneName]
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


def read_genes_asDictionary(genefile):
    '''Creates gene dictionary reader from gene file.'''
    return csv.DictReader(genefile, delimiter='\t')


def read_searchKeys(file):
    searchItems = []
    for line in file:
        term = line.strip()
        searchItems.append(term)
    return searchItems


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
        '--geneDb', metavar = 'GENEDB', type = file,
        help ='The tab separated file containing the gene information provided by ncbi, parsed representation of Hsxml')
    parser.add_argument(
        '--online', metavar = 'ONLINE', type = str, 
        help = 'yes: Search gene information from online (NCBI).')
    parser.add_argument(
        '--saveCache', metavar = 'CACHEFILE', type = str, 
        help = 'Save a cache of the downloaded results from NCBI into this file.')
    parser.add_argument(
        '--loadCache', metavar = 'LOADCACHE', type=str, 
        help = 'load previously cached download results from NCBI from this file.')
         
    return parser.parse_args()


def main():
    args = parse_args()
    genes = read_genes_asDictionary(args.genes)
    keywords = read_searchKeys(args.searchKeys)

    if args.geneDb != None:
        # search from the given db file
        geneDict = read_geneDb_asDictionary(args.geneDb)
        search_genes(genes, keywords, geneDict)
        
    elif args.online == 'yes':
        # get gene information over online
        fetchResults = fetch_geneInformation(genes)
        
        if args.saveCache != None:
            # save cache
            with open(args.saveCache, 'w') as file:
                pickle.dump(fetchResults, file)

        # searche keywords
        #search_genes_xml(genes, keywords, fetchResults)
    elif args.loadCache != None:
        # get gene information from cache
        try:
            cache = open(args.loadCache, 'r')
            fetchResults = pickle.load(cache)
            
            #search_genes_xml(genes, keywords, fetchResults)
        except:
            print 'Invalid cache file'

    else:  
        print 'Invalid argument'
        sys.exit()



  

if __name__ == '__main__':
    main()


