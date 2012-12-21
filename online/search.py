#!/bin/env python

import sys
import csv
from argparse import ArgumentParser

from Bio import Entrez
from lxml import etree

'''
Version modified from previous to accept Entrez_parser.py parsed xml file such as resulting from Homo_sapiens.xml. Accepts this output file, searchKeys and genes in tsv format. Works offline so that any number of variants can be tested.
'''


# for Entrez
Entrez.email = "sorik@student.unimelb.edu.au"

def get_geneIds(geneFileName):
    '''Gets gene Ids through ncbi search query.
       The genes are limited to Homo sapiens.
       'geneFileName' is a file including gene names to be searched. 
    '''
    geneIds = []
    geneFile = open(geneFileName, 'r')
    genes = csv.DictReader(geneFile, delimiter='\t')
    
    for gene in genes:
        searchTerm = 'Homo sapiens[organism] AND {}[sym]'.format(gene['Gene'])
        request = Entrez.esearch(db='gene', term=searchTerm)
        result = Entrez.read(request)
        geneIds.extend(result['IdList'])
    
    geneFile.close()
    return geneIds


def fetch_geneInformation(geneFileName):
    '''Searches gene information over online (NCBI).
       'genefile' is a file including gene names to be searched.
        1. get gene Ids 
        2. ePost using gene Ids -> get WebEnv and query_key
        3. eFetch using WebEnv and query_key -> get xml containing gene information
    '''
    
    geneIds = get_geneIds(geneFileName)

    idString = ','.join(geneIds)
    postRequest = Entrez.epost(db='gene', id=idString)
    postResult = Entrez.read(postRequest)
    
    webEnv = postResult['WebEnv']
    queryKey = postResult['QueryKey']
    fetchRequest = Entrez.efetch(db='gene',
                                 webenv=webEnv,
                                 query_key=queryKey,
                                 retmode='xml')
    # read in xml format
    fetchResults = fetchRequest.read()
    
    return fetchResults


# need to update
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

def merge_geneContent(geneContent, values):

    for key, value in values:
        geneContent[key] += value

    return geneContent


def get_geneName(geneEntry):
    name = geneEntry.xpath('./Entrezgene_gene/Gene-ref/Gene-ref_locus')
    if len(name) > 0:
        return name[0].text   
  
 
def get_otherSourceAnchor(entry):
    '''Find the element 'Other-source_anchor', and
       return the list of achor.
    '''
    result = []
    for child in entry.iterchildren():
        source = child.find('.//Gene-commentary_source')
        if source == None:
            continue
        for othersource in source.iterchildren():
            if othersource.tag == 'Other-source':
                anchor = othersource.find('.//Other-source_anchor')
                if anchor != None and anchor.text != None:
                    result.append(anchor.text)

    return result

def parse_geneCommentary(commentary):
    '''Parse element 'Gene-commentary'
    '''

    retValues = []
    # RIFs, pathway, PMID, interactions, conserved 
    rifs = []
    pathways = []
    pmids = []
    interactions = []
    conserved = []
    functions = []
    components = []
    processes = []

    type = ''
    heading = ''
    label = ''

    for item in commentary.iterchildren():        
        if item.tag == 'Gene-commentary_type':
            type = item.attrib['value']

        if item.tag == 'Gene-commentary_heading':
            heading = item.text

        if item.tag == 'Gene-commentary_label':
            label = item.text

        if item.tag == 'Gene-commentary_text':
            if type == 'generif' and item.text != None: 
                rifs.append(item.text)
                retValues.append(('GeneRIFs', rifs))
  
        if item.tag == 'Gene-commentary_refs':
            for pub in item.iterchildren():
                pubMedId = pub.find('.//Pub_pmid/PubMedId')
                if pubMedId != None and pubMedId != None:
                    pmids.append(pubMedId.text)
            retValues.append(('PmIds', pmids))

        if item.tag == 'Gene-commentary_products':
            for child in item.iterchildren():
                retValues += parse_geneCommentary(child)

        if item.tag == 'Gene-commentary_comment':
            if heading == 'Pathways':
                for child in item.iterchildren():
                    pathway = child.find('.//Gene-commentary_text')
                    if pathway != None and pathway.text != None:
                        pathways.append(pathway.text)
                retValues.append(('Pathways', pathways))

            elif heading == 'Interactions':
                interactions = []
                for child in item.iterchildren():
                    subComment = item.find('.//Gene-commentary_comment')
                    for grandChild in child.iterchildren():
                        result = get_otherSourceAnchor(grandChild)
                        interactions += result
                retValues.append(('Interactions', interactions))
                 
            elif heading == 'NCBI Reference Sequences (RefSeq)':
                for child in item.iterchildren():
                    retValues += parse_geneCommentary(child)

            elif heading == 'Conserved Domains':
                result = get_otherSourceAnchor(item)
                retValues.append(('Conserved Domains', result))  
                 
            elif label == 'Function':
                result = get_otherSourceAnchor(item)
                retValues.append(('Function', result))    

            elif label == 'Component':
                result = get_otherSourceAnchor(item)
                retValues.append(('Component', result))  

            elif label == 'Process':
                result = get_otherSourceAnchor(item)
                retValues.append(('Process', result))  

            else:
                #Finding 'Conserved Domains'
                for child in item.iterchildren():
                    if child.tag == 'Gene-commentary':
                        heading = child.find('.//Gene-commentary_heading')
                        if heading != None and heading.text == 'Conserved Domains':
                            retValues += parse_geneCommentary(child)
                 
    return retValues
    


def get_geneContent(geneEntry):
    '''Parse element Entrezgene. 
    '''

    geneContent = {'Gene Name' : [], 
                   'Description' : [],
                   'Synonyms' : [],
                   'Alternative Name' : [],
                   'Summary' : [],
                   'GeneRIFs' : [],
                   'PmIds' : [],
                   'Pathways' : [],
                   'Interactions' : [],
                   'Conserved Domains' : [],
                   'Function' : [],
                   'Component' : [],
                   'Process' : []
                  }

    for elem in geneEntry.iterchildren():

        if elem.tag == 'Entrezgene_track-info':
            status = elem.find('.//Gene-track/Gene-track_status') 
            if status != None and status.attrib['value'] == 'discontinued':
                continue

        elif elem.tag == 'Entrezgene_gene':
            ref = elem.find('.//Gene-ref')
            name = ref.find('.//Gene-ref_locus')
            desc = ref.find('.//Gene-ref_desc')
            syn = ref.find('.//Gene-ref_syn')
            if name != None and name.text != None: 
                geneContent['Gene Name'].append(name.text)
            if desc != None and desc.text != None: 
                geneContent['Description'].append(desc.text)
            if syn != None:
                for child in syn.iterchildren():
                    if child.tag == 'Gene-ref_syn_E' and child.text != None:
                        geneContent['Synonyms'].append(child.text)   

        elif elem.tag == 'Entrezgene_summary' and elem.text != None:
            geneContent['Summary'].append(elem.text)
            
        elif elem.tag == 'Entrezgene_comments':
            for subElem in elem.iterchildren():
                result = parse_geneCommentary(subElem)
                geneContent = merge_geneContent(geneContent, result)

        elif elem.tag == 'Entrezgene_properties':
            entry = elem.find('.//Gene-commentary/Gene-commentary_comment')
            if entry != None:
                for child in entry.iterchildren():
                    result = parse_geneCommentary(child)
                    geneContent = merge_geneContent(geneContent, result)

        elif elem.tag == 'Entrezgene_prot':
            alterName = elem.find('.//Prot-ref/Prot-ref_name')
            if alterName != None:
                for name in alterName.iterchildren():
                    if name.tag == 'Prot-ref_name_E' and name.text != None:
                        geneContent['Alternative Name'].append(name.text)

    
    return geneContent
            

def search_genes_xml(genefilename, keywords, topNRank, xmlfile):
    '''While parsing xml file, search keywords from the xml contents.
       If the keyword is found, it adds the additional column to the genefile.
       The additional column describes the keyword, its rank, and where the keyword is found.
    '''
     
    genes_searched = {}
    additionalFields = ['Keyword, Rank, DB Fields']

    # read genes and extract gene names which we want to find 
    geneNames = []
    genefile = open(genefilename, 'r')
    genes = csv.DictReader(genefile, delimiter='\t')
    for gene in genes: 
        if gene['Func'] != 'intergenic':
            geneNames.append(gene['Gene'])
    genefile.close()

    try:
        content = etree.iterparse(xmlfile, events=('end',), tag='Entrezgene')
    except: # in the case when reading xml is failed
        print 'Fail to read xml file (%s).'%xmlfile
        sys.exit()

    for event, geneEntry in content:
        # get gene name of geneEntry, and 
        # if the gene name is on the list of gene names which are in gene file,
        # parse the gene content from xml and search keywords over the parsed content
        geneName = get_geneName(geneEntry) 
        
        if geneName in geneNames: 
            dbEntry = get_geneContent(geneEntry)
            searchResult = search_keywords_inDict(dbEntry, keywords, topNRank)

            resultStr = ''
            for keyword, rank, fields  in searchResult:
                resultStr += '('+','.join([keyword, str(rank), ';'.join(fields)])+')'
            genes_searched[geneName] = resultStr

        # clear node references
        geneEntry.clear()
        while geneEntry.getprevious() is not None:
           del geneEntry.getparent()[0]

    del content

    # print the results
    genefile = open(genefilename, 'r')
    genes = csv.DictReader(genefile, delimiter='\t')    
    headers = genes.fieldnames 
    newHeaders = headers + additionalFields
    print '\t'.join(newHeaders)

    for gene in genes:
        geneInfo = [gene[field] for field in headers]
        geneName = gene['Gene']
        try:
            geneInfo.append(genes_searched[geneName])        
        except KeyError:   # in case when keyword is not found
            geneInfo.append('')
        print '\t'.join(geneInfo)

    genefile.close()


      
def search_genes(genefilename, keywords, geneDict):
    '''Searches gene information from geneDict, and
       prints the results.
           genes: csv.DictReader, keywords: list, geneDict: dict
       1. searches keywords, and appends results to the original data
          ['Keyword', 'Key Rank']
    '''
    additionalFields = ['Keyword', 'Key Rank']
   
    genefile = open(genefilename, 'r')
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
                dbEntry = geneDict[geneName]
                searchResults = search_keywords_inString(dbEntry, keywords)
                if searchResults != None:
                    keyword, rank = searchResults
                    geneInfo += [keyword, str(rank)]
                else:
                    # if db does not contain the keywords
                    geneInfo += ['', '']
            except KeyError:
                # if db does not contain the geneName
                geneInfo += ['', '']
     
        print '\t'.join(geneInfo)

    genefile.close()

   
def search_keywords_inDict(dbEntry, keywords, topNRank):
    '''Search top N ranking keywords from dbEntry.
       dbEntry is a dictionary which contains information of each field.
       The value of each key is a list.
       e.g) {'AlterName' : ['abc', 'def'], ...}
       Returns a list of tuples containging 
                         the keyword hitted, 
                         the rank of keyword, and 
                         the fields where the keyword hitted.
    '''
    keysFound = []
    fields = []
    numFound = 0
    for n, keyword in enumerate(keywords):
        for key in dbEntry:
            for content in dbEntry[key]:
                if keyword in content:
                    if key not in fields:
                        fields.append(key)
        if len(fields) > 0:
            keysFound.append((keyword, n+1, fields))
            fields = []
            numFound += 1
        if numFound >= topNRank:
            break
   
    return keysFound
         

def search_keywords_inString(dbEntry, keywords):
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
        '--topNRank', metavar = 'TOPNRANK', type = int, default=1,
        help = 'Top N rank would be searched. Default is 1.') 
    parser.add_argument(
        '--genes', metavar = 'GENES', type = str, required=True,
        help ='The tab separated file containing the gene information including name of the gene, one gene name per line')
    parser.add_argument(
        '--geneDb', metavar = 'GENEDB', type = file,
        help ='The tab separated file containing the gene information provided by ncbi, parsed representation of Hsxml')
    parser.add_argument(
        '--online', action='store_true', 
        help = 'Search gene information from online (NCBI).')
    parser.add_argument(
        '--saveCache', metavar = 'CACHEFILE', type = str, 
        help = 'Save a cache of the downloaded results from NCBI into this file.')
    parser.add_argument(
        '--loadCache', metavar = 'LOADCACHE', type=str, 
        help = 'load previously cached download results from NCBI from this file.')
         
    return parser.parse_args()


def main():
    args = parse_args()
    keywords = read_searchKeys(args.searchKeys)

    if args.geneDb != None:
        # search from the given db file
        geneDict = read_geneDb_asDictionary(args.geneDb)
        search_genes(args.genes, keywords, geneDict)
        
    elif args.online == True: 
        # get gene information over online
        fetchResults = fetch_geneInformation(args.genes)
        
        cacheFile = ''
        if args.saveCache != None:
            cacheFile = args.saveCache
        else:
            cacheFile = 'cache'
        # save cache
        with open(cacheFile, 'w') as file:
            file.write(fetchResults)
            file.close()
            
        # searche keywords
        search_genes_xml(args.genes, keywords, args.topNRank, cacheFile)
    
    elif args.loadCache != None:
        # get gene information from cache
        search_genes_xml(args.genes, keywords, args.topNRank, args.loadCache)

    else:  
        print 'Invalid argument'
        sys.exit()

  
   
    

  

if __name__ == '__main__':
    main()


