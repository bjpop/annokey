#!/bin/env python

'''
--------------------------------------------------------------------------------

NCBI Gene Search Tool
---------------------

Authors:   Daniel Park, Sori Kang, Bernie Pope, Tu Nguyen-Dumont.
Copyright: 2013
Website:   https://github.com/bjpop/annokey
License:   BSD, see LICENCE file in source distribution. 

This program searches NCBI database for genes of interest
based on concept-keyword search.

This program supports both online search and offline search.
If the option --online is on, this program downloads gene database in xml
format from the NCBI server and search keywords. The option --saveCache allows
the downloaded information to be saved int the user's directory.

For offline search, there are two options. --loadCache and --geneDb.
--loadCache accepts a xml file for gene information, while --geneDb accepts
a tab separated file that containing gene information.

The search results are appended at the end of column of the
input gene file with the information about the keyword hit,
the keyword rank, and the section where the keyword hit.

Required inputs:

    --searchKeys FILENAME   a text file of search keywords/keyphrases. One
                            term per line.

    --genes FILENAME        a text file of gene information, one gene per line.
                            Format is tab separated. Must contain at least one
                            column of gene names.

--------------------------------------------------------------------------------
'''

import sys
import csv
import os
import glob
from argparse import ArgumentParser
from Bio import Entrez
from lxml import etree
from StringIO import StringIO
import itertools
import cPickle as pickle

# for Entrez
# XXX Have to consider registering tool and email.
# XXX this should be a parameter.
Entrez.email = "bjpop@.unimelb.edu.au"

#NCBI access rate limit
#http://www.ncbi.nlm.nih.gov/books/NBK25497/
#In order not to overload the E-utility servers, NCBI recommends that users
#post no more than three URL requests per second and limit large jobs to either
#weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays. Failure
#to comply with this policy may result in an IP address being blocked from
#accessing NCBI. If NCBI blocks an IP address, service will not be restored
#unless the developers of the software accessing the E-utilities register values
#of the tool and email parameters with NCBI. The value of tool should be a
#string with no internal spaces that uniquely identifies the software producing
#the request. The value of email should be a complete and valid e-mail address
#of the software developer and not that of a third-party end user. The value of
#email will be used only to contact developers if NCBI observes requests that
#violate our policies, and we will attempt such contact prior to blocking
#access. In addition, developers may request that the value of email be added
#to the E-utility mailing list that provides announcements of software updates,
#known bugs and other policy changes affecting the E-utilities. To register tool
#and email values, simply send an e-mail to eutilities@ncbi.nlm.nih.gov
#including the desired values along with the name of either a developer or the
#organization creating the software. Once NCBI establishes communication with a
#developer, receives values for tool and email and validates the e-mail address
#in email, the block will be lifted. Once tool and email values are registered,
#all subsequent E-utility requests from that software package should contain
#both values. Please be aware that merely providing values for tool and email
#in requests is not sufficient to comply with this policy; these values must be
#registered with NCBI. Requests from any IP that lack registered values for tool
#and email and that violate the above usage policies may be blocked. Software
#developers may register values of tool and email at any time, and are
#encouraged to do so.

# The column positions of gene name and func.
# XXX can we get rid of these?
gene_name_pos = None
func_pos = None
has_column_headers = None

def set_column_options(args):
    '''Find the column index containing gene name or func.
       If both xxCol and xxColPos are presented, then use xxxCol.
    '''
    global gene_name_pos, gene_func_pos, has_column_headers
    # We assume that if xxColPos is presented,
    # the file does not have a column headers.
    # XXX is this assumption reasonable?
    if args.geneNameColPos and args.geneFuncColPos:
        has_column_headers = False
        gene_name_pos = args.geneNameColPos
        gene_func_pos = args.geneFuncColPos
    else:
        has_column_headers = True
        with open(args.genes, 'rb') as file:
            genefile = csv.reader(file, delimiter='\t')
            column_headers = genefile.next()
            name_column = args.geneNameCol if args.geneNameCol else 'Gene'
            func_column = args.geneFuncCol if args.geneFuncCol else 'Func'
            try:
                if args.geneNameColPos:
                    gene_name_pos = args.geneNameColPos
                else:
                    gene_name_pos = column_headers.index(name_column)
            except ValueError:
                print '* Error * '
                print 'Could not find column for gene name (%s).' % name_column
                print 'Column header for gene name should be Gene.'
                print 'Or use --geneNameCol or --geneNameColPos options.'
                return False

            try:
                if args.geneFuncColPos:
                    gene_func_pos = args.geneFuncColPos
                else:
                    gene_func_pos = column_headers.index(func_column)
            except ValueError:
                print '* Error *'
                print 'Could not find column for func (%s).' % func_column
                print 'Column header for gene func should be Func.'
                print 'Or use --geneFuncCol or --geneFuncColPos options.'
                return False
    return True


def get_column_headers(genefileName):
    '''Read column headers if the file has it'''
    if has_column_headers:
        with open(genefileName, 'rb') as file:
            genefile = csv.reader(file, delimiter='\t')
            return genefile.next()


def read_gene_file(genefileName, padding=False):
    '''Read gene file and return generator.'''

    max_column_length = 0
    if padding:
        # Find the maximum column length
        with open(genefileName, 'rb') as file:
            genefile = csv.reader(file, delimiter='\t')
            for line in genefile:
                if len(line) > max_column_length:
                    max_column_length = len(line)

    with open(genefileName, 'rb') as file:
        genefile = csv.reader(file, delimiter='\t')
        if has_column_headers:
            genefile.next()
        for line in genefile:
            # gene name including (,;
            # If '(' included, then use the gene name before '('
            # XXX If ';' included, then use the first gene name.
            # Need to consider use both(first and second)  gene names.
            # http://www.openbioinformatics.org/annovar/annovar_gene.html
            func = line[gene_func_pos].strip()
            name = line[gene_name_pos].strip()
            if func == 'splicing':
                pa_pos = name.find('(')
                if pa_pos != -1:
                    name = name.split('(')[0]
            if func == 'exonic;splicing':
                se_pos = name.find(';')
                if se_pos != -1:
                    name = name.split(';')[0]
            if padding:
                num_to_add = max_column_length - len(line)
                line.extend(['' for i in range(num_to_add)])
            yield name, func, line


def get_geneIds(genefilename):
    '''Get gene Ids through NCBI search query.'''

    def get_gene_names(genefilename, count):
        names = []
        for n, (name, func, info) in enumerate(read_gene_file(genefilename)):
            names.append(name)
            if (n+1) % count == 0:
                yield names
                names = []
        yield names

    # Send request to server for every 100 genes.
    geneIds = []
    for names in get_gene_names(genefilename, 100):
        terms = ' OR '.join(['%s[sym]' % name for name in names])
        search_term = 'Homo sapiens[organism] AND (%s)' % terms
        request = Entrez.esearch(db='gene',
                                 term=search_term,
                                 retmax=10000)
        result = Entrez.read(request)
        geneIds.extend(result['IdList'])

    geneIds = make_unique_list(geneIds)
    return geneIds


def make_unique_list(list):
    item_seen = {}
    for i in list:
        item_seen[i] = 1
    return item_seen.keys()


def lookup_pubmed_ids(ids):
    not_cached_ids = []
    # search for all the cached pubmed ids first, and
    # collect the non-cached ones.
    for id in ids:
        cache_result = lookup_pubmed_cache(id)
        if cache_result is not None:
            #print("found {} in cache:".format(id))
            yield cache_result
        else:
            print("did not find {} in cache:".format(id))
            not_cached_ids.append(id)

    # I don't think we should do the chunking here, but instead
    # rely on the Entrz history.
    # fetch the non-cached ids from NCBI

    if len(not_cached_ids) == 0:
        return

    idString = ','.join(not_cached_ids)
    postRequest = Entrez.epost(db='pubmed', id=idString)
    postResult = Entrez.read(postRequest)
    webEnv = postResult['WebEnv']
    queryKey = postResult['QueryKey']
    retstart = 0
    chunk_size = 10000

    while retstart < len(not_cached_ids):
        try:
            fetch_request = Entrez.efetch(db='pubmed',
                                          webenv=webEnv,
                                          query_key=queryKey,
                                          retmode='xml',
                                          retmax=chunk_size,
                                          retstart=retstart)
            pubmed_result = fetch_request.read()
        # XXX What should we do on exception? Try again? Sleep? Quit?
        except Exception as e:
            print("fetch failed with: {} {}".format(e, type(e)))
            break 

        content = etree.iterparse(StringIO(pubmed_result), events=('end',), tag='PubmedArticle')
        for event, pubmed_entry in content:
            # XXX we should cache this result possibly
            yield get_pubmed_content(pubmed_entry)
            # Clear node references
            pubmed_entry.clear()
            while pubmed_entry.getprevious() is not None:
                del pubmed_entry.getparent()[0]
        del content
        retstart += chunk_size


def lookup_pubmed_cache(id):
    hashed_id = int(id) % 256
    pickle_filename = os.path.join('pubmed_cache', str(hashed_id), id)
    if os.path.isfile(pickle_filename):
        with open(pickle_filename) as file:
            pubmed_content = pickle.load(file)
            return pubmed_content
    else:
        return None

def fetch_and_save(db, ids, filename):
    '''Fetch using Entrez'''
    # Before posting, make a list of unique Ids
    uniqIds = make_unique_list(ids)
    idString = ','.join(uniqIds)
    postRequest = Entrez.epost(db=db, id=idString)
    postResult = Entrez.read(postRequest)

    webEnv = postResult['WebEnv']
    queryKey = postResult['QueryKey']
    fetchRequest = Entrez.efetch(db=db,
                                 webenv=webEnv,
                                 query_key=queryKey,
                                 retmode='xml')
    fetchResults = fetchRequest.read()
    try:
        with open(filename, 'w') as file:
            file.write(fetchResults)
    except (KeyboardInterrupt, IOError):
        if os.path.exists(filename):
            os.remove(filename)
            raise
    return filename

# What to be saved is
# <PubmedArticle
#              <MedlineCitation
#                              <PMID
#                              <Article
#                                      <Journal
#                                              <JournalIssue
#                                                           <Volume
#                                                           <Issue
#                                                           <PubDate
#                                                                   <Year
#                                              <Title
#                                     <ArticleTitle
#                                     <Pagination
#                                                <MedlinePgn
#                                     <Abstract
#                                              <AbstractText
#                                               Label='BACKGROUND'
#                                              <AbstractText
#                                               Label='MATERIAL AND METHODS'
#                                              <AbstractText
#                                               Label='RESULTS'
#                                              <AbstractText
#                                               Label='CONCLUSION'
#


def get_pubmed_content(pubmedEntry):
    '''Parse element PubmedArticle'''
    pubmedContent = {}
    for elem in pubmedEntry.iterchildren():
        if elem.tag == 'MedlineCitation':
            pmid = elem.find('.//PMID')
            if pmid is not None and pmid.text is not None:
                pubmedContent['PMID'] = pmid.text
            article = elem.find('.//Article')
            if article is not None:
                journal = elem.find('.//Journal')
                articleTitle = elem.find('.//ArticleTitle')
                pagination = elem.find('.//Pagination')
                abstract = elem.find('.//Abstract')
            if journal is not None:
                volume = journal.find('.//JournalIssue/Volume')
                issue = journal.find('.//JournalIssue/Issue')
                pubdate = journal.find('.//JournalIssue/PubDate/Year')
                journalTitle = journal.find('.//Title')
            if volume is not None and volume.text is not None:
                pubmedContent['Volume'] = volume.text
            if issue is not None and issue.text is not None:
                pubmedContent['Issue'] = issue.text
            if pubdate is not None and pubdate.text is not None:
                pubmedContent['Year'] = pubdate.text
            if journalTitle is not None and journalTitle.text is not None:
                pubmedContent['Journal Title'] = journalTitle.text
            if articleTitle is not None and articleTitle.text is not None:
                pubmedContent['Article Title'] = articleTitle.text
            if pagination is not None and pagination.text is not None:
                pgn = pagination.find('.//MedlinePgn')
                if pgn is not None and pgn.text is not None:
                    pubmedContent['Pagination'] = pgn.text
            if abstract is not None:
                abstractText = ''
                for abstractChild in abstract.iterchildren():
                    if abstractChild.tag == 'AbstractText':
                        if abstractChild.text is not None:
                            abstractText += abstractChild.text
                pubmedContent['Abstract'] = abstractText
    return pubmedContent

# XXX need to update
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
    '''Merge gene information.

    Args:
         geneContent: dict containing gene information.
                      values to be merged into this geneContent.
         values: a list of tuple containing gene information.
                 It is to be merged into geneContent.
    Returns:
         A dict containing gene information.
    '''

    for key, value in values:
        geneContent[key] += value 
    return geneContent


def get_geneName(geneEntry):
    name = geneEntry.find('.//Entrezgene_gene/Gene-ref/Gene-ref_locus')
    if name is not None and name.text is not None:
        return name.text


def get_otherSourceAnchor(entry):
    '''Find the element 'Other-source_anchor',
       and return the list of achor.
    '''
    result = []
    for child in entry.iterchildren():
        source = child.find('.//Gene-commentary_source')
        if source is not None:
            for othersource in source.iterchildren():
                if othersource.tag == 'Other-source':
                    anchor = othersource.find('.//Other-source_anchor')
                    if anchor is not None and anchor.text is not None:
                        result.append(anchor.text)

    return result


def parse_geneCommentary(commentary):
    '''Parse element 'Gene-commentary'''

    retValues = []

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
            if type == 'generif' and item.text is not None:
                rifs.append(item.text)
                retValues.append(('GeneRIFs', rifs))

        if item.tag == 'Gene-commentary_refs':
            for pub in item.iterchildren():
                pubMedId = pub.find('.//Pub_pmid/PubMedId')
                if pubMedId is not None and pubMedId.text is not None:
                    pmids.append(pubMedId.text)
            retValues.append(('PmIds', pmids))

        if item.tag == 'Gene-commentary_products':
            for child in item.iterchildren():
                retValues += parse_geneCommentary(child)

        if item.tag == 'Gene-commentary_comment':
            if heading == 'Pathways':
                for child in item.iterchildren():
                    pathway = child.find('.//Gene-commentary_text')
                    if pathway is not None and pathway.text is not None:
                        pathways.append(pathway.text)
                retValues.append(('Pathways', pathways))

            elif heading == 'Interactions':
                interactions = []
                for child in item.iterchildren():
                    subComment = item.find('.//Gene-commentary_comment')
                    for grandChild in child.iterchildren():
                        interactions += get_otherSourceAnchor(grandChild)
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
                # Find 'Conserved Domains'
                for child in item.iterchildren():
                    if child.tag == 'Gene-commentary':
                        heading = child.find('.//Gene-commentary_heading')
                        if (heading is not None and
                                heading.text == 'Conserved Domains'):
                            retValues += parse_geneCommentary(child)

    return retValues


def get_geneContent(geneEntry):
    '''Parse element Entrezgene.'''

    geneContent = {'Gene Name': [],
                   'Description': [],
                   'Synonyms': [],
                   'Alternative Name': [],
                   'Summary': [],
                   'GeneRIFs': [],
                   'PmIds': [],
                   'Pathways': [],
                   'Interactions': [],
                   'Conserved Domains': [],
                   'Function': [],
                   'Component': [],
                   'Process': []}

    for elem in geneEntry.iterchildren():

        if elem.tag == 'Entrezgene_track-info':
            status = elem.find('.//Gene-track/Gene-track_status')
            if (status is not None and
                    status.attrib['value'] == 'discontinued'):
                continue

        elif elem.tag == 'Entrezgene_gene':
            ref = elem.find('.//Gene-ref')
            name = ref.find('.//Gene-ref_locus')
            desc = ref.find('.//Gene-ref_desc')
            syn = ref.find('.//Gene-ref_syn')
            if name is not None and name.text is not None:
                geneContent['Gene Name'].append(name.text)
            if desc is not None and desc.text is not None:
                geneContent['Description'].append(desc.text)
            if syn is not None:
                for child in syn.iterchildren():
                    if (child.tag == 'Gene-ref_syn_E' and
                            child.text is not None):
                        geneContent['Synonyms'].append(child.text)

        elif elem.tag == 'Entrezgene_summary' and elem.text is not None:
            geneContent['Summary'].append(elem.text)

        elif elem.tag == 'Entrezgene_comments':
            for subElem in elem.iterchildren():
                result = parse_geneCommentary(subElem)
                geneContent = merge_geneContent(geneContent, result)

        elif elem.tag == 'Entrezgene_properties':
            entry = elem.find('.//Gene-commentary/Gene-commentary_comment')
            if entry is not None:
                for child in entry.iterchildren():
                    result = parse_geneCommentary(child)
                    geneContent = merge_geneContent(geneContent, result)

        elif elem.tag == 'Entrezgene_prot':
            alterName = elem.find('.//Prot-ref/Prot-ref_name')
            if alterName is not None:
                for name in alterName.iterchildren():
                    if (name.tag == 'Prot-ref_name_E' and
                            name.text is not None):
                        geneContent['Alternative Name'].append(name.text)

    return geneContent


def search_genes_in_xml(genefilename, keywords, topNRank, genexmlfile):
    '''While parsing xml file, search keywords from the xml contents.
       If the keyword is found, it adds the additional column to the genefile.
       The additional column describes
       the keyword, its rank, and where the keyword is found.
    '''
    def get_keyword_hits_from_pubmed(pmIds):
        '''Return pubmed_hits'''
        pubmed_hit = {}
        for pubmed_content in lookup_pubmed_ids(pmIds):
            #print("pubmed for {}".format(pubmed_content['PMID']))
            #print(pubmed_content)
            keywords_found = search_keywords_inPubMed(pubmed_content, keywords, topNRank)
            pubmed_hit[pubmed_content['PMID']] = keywords_found
        return pubmed_hit

    def get_keyword_hits_from_genexml(genefilename, genexmlfile):
        '''Return gene_searched and total PMIDs'''
        genes_searched = {}
        pmIds = set()
        geneNames = []
        # Read genes and extract gene names which we want to find.
        for name, func, info in read_gene_file(genefilename):
            if func != 'intergenic':
                geneNames.append(name)

        content = etree.iterparse(genexmlfile, events=('end',), tag='Entrezgene')
        for event, geneEntry in content:
            # Get gene name of geneEntry, and
            # if the gene name is on the list of gene names which
            # are in gene file, parse the gene content from xml and
            # search keywords over the parsed content.
            geneName = get_geneName(geneEntry)
            if geneName in geneNames:
                dbEntry = get_geneContent(geneEntry)
                searchResult = search_keywords_inDict(dbEntry, keywords, topNRank)
                if geneName in genes_searched:
                    record = genes_searched[geneName]
                    record['keywordHit'] += searchResult
                    # XXX this should be a set
                    record['PmIds'] += make_unique_list(dbEntry['PmIds'])
                else:
                    results = {}
                    results['keywordHit'] = searchResult
                    results['PmIds'] = make_unique_list(dbEntry['PmIds'])
                    genes_searched[geneName] = results
                pmIds.update(results['PmIds'])
            # Clear node references
            geneEntry.clear()
            while geneEntry.getprevious() is not None:
                del geneEntry.getparent()[0]
        del content
        return genes_searched, pmIds


    def merge_gene_and_pubmed_hits(genes_searched, pubmed_hit, topNRank):
        '''Return genes_searched including pubmed_hits'''
        # Match gene and pubmed
        for gene, info_searched in genes_searched.iteritems():
            gene_pubmed_hit = {}
            gene_keyword_hit = info_searched['keywordHit']
            gene_pmids = info_searched['PmIds']
            # Find pubmed hit information over all PMIDs related to the gene.
            # Count the number of hits of the keyword over all PMIDs.
            for pmid in gene_pmids:
                try:
                    for keyword, rank in pubmed_hit[pmid]:
                        if keyword in gene_pubmed_hit:
                            rank, count = gene_pubmed_hit[keyword]
                            gene_pubmed_hit[keyword] = (rank, count+1)
                        else:
                            gene_pubmed_hit[keyword] = (rank, 1)
                except KeyError:
                    pass
            # Merge keyword sections if the keyword was also found in PMID.
            for n, (keyword, rank, fields) in enumerate(gene_keyword_hit):
                if keyword in gene_pubmed_hit:
                    rank, count = gene_pubmed_hit[keyword]
                    fields.append('PMID(%s/%s)' % (count, len(gene_pmids)))
                    gene_keyword_hit[n] = (keyword, rank, fields)
                    del gene_pubmed_hit[keyword]
            # For the remaining PMID hits, append to keywordHit.
            for keyword, (rank, count) in gene_pubmed_hit.items():
                fields = ['PMID(%s/%s)' % (count, len(gene_pmids))]
                gene_keyword_hit.append((keyword, rank, fields))
            # Sort keywordHit according to its rank,
            # so we can determine top n rank.
            gene_keyword_hit.sort(key=lambda tup: tup[1])
            genes_searched[gene]['keywordHit'] = gene_keyword_hit[:topNRank]
        return genes_searched

    genes_searched, pmIds = get_keyword_hits_from_genexml(genefilename,
                                                          genexmlfile)

    pubmed_hit = get_keyword_hits_from_pubmed(pmIds)

    genes_searched = merge_gene_and_pubmed_hits(genes_searched,
                                                pubmed_hit, topNRank)


    # Make a string for printing.
    # If there is no keyword hit, then the output would be ''
    # Otherwise, each keyword would be concatenated with ().
    for name, info_searched in genes_searched.iteritems():
        output_format = ''
        for keyword, rank, fields in info_searched['keywordHit']:
            output_format += '(%s,%s,%s)' % (keyword, rank, ';'.join(fields))
        genes_searched[name]['keywordHit'] = output_format

    # Print the results.
    additionalFields = ['Keyword,Rank,Sections']
    headers = get_column_headers(genefilename)
    if headers:
        headers += additionalFields
        print '\t'.join(headers)
    for name, func, geneInfo in read_gene_file(genefilename, padding=True):
        if name in genes_searched:
            annotation = genes_searched[name]['keywordHit']
        else:
            annotation = "Could not find '%s' in NCBI database" % name
        if annotation != '':
            print '\t'.join(geneInfo + [annotation])


def search_genes(genefilename, keywords, topNRank, geneDict):
    '''Search gene information from geneDict, and
       print the results.
           genes: csv.DictReader, keywords: list, geneDict: dict
       1. searches keywords, and appends results to the original data
          ['Keyword,Rank']
    '''
    column_headers = get_column_headers(genefilename)
    if column_headers:
        column_headers.append('Keyword,Rank')
        print '\t'.join(column_headers)

    for name, func, geneInfo in read_gene_file(genefilename, padding=True):
        if func == 'intergenic':
            geneInfo.append('')
            continue
        # Search keywords against gene database.
        result = ''
        try:
            dbEntry = geneDict[name]
            searchResults = search_keywords_inString(dbEntry,
                                                     keywords, topNRank)
            if searchResults:
                for key, rank in searchResults:
                    result += '(%s,%s)' % (key, rank)
            else:
                # XXX need detail information?
                result = ''
        except KeyError:
            # If DB does not contain the information on gene name.
            result = "Could not find '%s' in DB" % name
            print(geneDict.keys())
            exit()
        geneInfo.append(result)
        print '\t'.join(geneInfo)


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
    for n, keyword in enumerate(keywords):
        for field, content in dbEntry.iteritems():
            if field == 'GeneRIFs':
                hit = [keyword in item for item in content]
                if sum(hit) > 0:
                    fields.append('%s(%s/%s)' %
                                  (field, sum(hit), len(content)))
            else:
                for item in content:
                    if keyword in item and field not in fields:
                        fields.append(field)
        if len(fields) > 0:
            keysFound.append((keyword, n+1, fields))
            fields = []

        if len(keysFound) >= topNRank:
            break

    return keysFound


def search_keywords_inPubMed(pubmedContent, keywords, topNRank):
    '''Search keywords from pubmedContent.'''
    # Only interested in article title and abstract
    content = ''
    try:
        content += pubmedContent['Article Title']
        content += ';' + pubmedContent['Abstract']
    except KeyError:
        pass

    return search_keywords_inString(content, keywords, topNRank)


def search_keywords_inString(dbEntry, keywords, topNRank):
    '''Searche keyword in the order (according to the rank) and
       return the keyword and its rank (position in the list).
       If fail to search, (python) returns None.
       dbEntry is a string.
    '''
    keysFound = []
    for n, item in enumerate(keywords):
        if item in dbEntry:
            keysFound.append((item, n+1))
        if len(keysFound) >= topNRank:
            break
    return keysFound


def read_geneDb_asDictionary(geneDb):
    '''Create gene dictionary from gene database.
       geneDb is a string in which the content of genes are separated by tab.
    '''
    geneDict = {}
    for line in geneDb:
        geneContent = line.split('\t')
        # .tsv file has no column headers.
        # index 1(one) is gene name.
        # there is a case where a line is empty.
        if len(geneContent) >= 2:
            geneDict[geneContent[1]] = line

    return geneDict


def parse_args():
    parser = ArgumentParser(description='Search NCBI for genes of interest, '
                                        'based on concept-keyword search.')
    parser.add_argument('--searchKeys',
                        metavar='SEARCHKEYS',
                        type=file,
                        required=True,
                        help='The tab separated file containing '
                             'the keywords to be searched.')
    parser.add_argument('--topNRank',
                        metavar='TOPNRANK',
                        type=int,
                        default=1,
                        help='Top N rank that would be searched. '
                             'Default is 1.')
    parser.add_argument('--genes',
                        metavar='GENES',
                        type=str,
                        required=True,
                        help='The tab separated file containing '
                             'the gene information including '
                             'name of the gene, one gene name per line.')
    parser.add_argument('--geneDb',
                        metavar='GENEDB',
                        type=file,
                        help='The tab separated file containing '
                             'the gene information provided by ncbi, '
                             'parsed representation of Homo sapiens '
                             'gene information.')
    parser.add_argument('--online',
                        action='store_true',
                        help='Search gene information from online (NCBI).')
    parser.add_argument('--saveCache',
                        metavar='CACHEFILE',
                        type=str,
                        help='Save a cache of the downloaded results '
                             'from NCBI into this file.')
    parser.add_argument('--loadCache',
                        metavar='CACHEFILE',
                        type=str,
                        help='Load previously cached, downloaded results '
                             'from NCBI.')
    parser.add_argument('--geneNameCol',
                        metavar='NANE',
                        type=str,
                        help='The name of the column containing gene name.'
                             'If the input gene file has column names and '
                             'if the column name indicating  gene name '
                             "is 'Gene', then there is no need to pass "
                             "this option. This program would use 'Gene' "
                             'as default to find the column that has '
                             'gene names')
    parser.add_argument('--geneNameColPos',
                        metavar='POSITION',
                        type=int,
                        help='The index of the column containing gene name. '
                             'The index should start from 0.'
                             'If the input gene file does not have column '
                             'names, then you should pass this option to '
                             'specify which column should be read for gene '
                             'names')
    parser.add_argument('--geneFuncCol',
                        metavar='NANE',
                        type=str,
                        help='The name of the column containing gene func.'
                             'If the input gene file has column names and '
                             'if the column name indicating  gene name '
                             "is 'Func', then there is no need to pass "
                             "this option. This program would use 'Func' "
                             'as default to find the column that has '
                             'gene names')
    parser.add_argument('--geneFuncColPos',
                        metavar='POSITION',
                        type=int,
                        help='The index of the column containing gene func. '
                             'The index should start from 0.'
                             'If the input gene file does not have column '
                             'names, then you should pass this option to '
                             'specify which column should be read for gene '
                             'func')
    return parser.parse_args()


def main():
    args = parse_args()
    # Set the column positions of gene name and func.
    # If the input arguments are not correct,
    # prints error message and finish the program.
    result = set_column_options(args)
    if not result:
        return

    keywords = [line.strip() for line in args.searchKeys]
    if args.geneDb:
        # Search from the given DB file.
        geneDict = read_geneDb_asDictionary(args.geneDb)
        search_genes(args.genes, keywords, args.topNRank, geneDict)

    elif args.online:
        cacheFile = ''
        if args.saveCache:
            cacheFile = args.saveCache
        else:
            cacheFile = 'cache'
        # Get gene information over online.
        geneIds = get_geneIds(args.genes)
        filename = fetch_and_save('gene', geneIds, cacheFile)
        # Search keywords.
        search_genes_in_xml(args.genes, keywords,
                            args.topNRank, filename)
    elif args.loadCache:
        # Get gene information from cache.
        search_genes_in_xml(args.genes, keywords,
                            args.topNRank, args.loadCache)

    else:
        print '* Error *'
        print 'One of --geneDb, --online, and --loadCache is necessary.'

if __name__ == '__main__':
    main()
