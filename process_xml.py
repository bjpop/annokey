#!/bin/env python

'''
--------------------------------------------------------------------------------

Annokey: a NCBI Gene Database Keyword Search Tool
-------------------------------------------------

Authors:   Daniel Park, Sori Kang, Bernie Pope, Tu Nguyen-Dumont.
Copyright: 2013
Website:   https://github.com/bjpop/annokey
License:   BSD, see LICENCE file in source distribution. 


--------------------------------------------------------------------------------
'''

from Bio import Entrez
from lxml import etree
from StringIO import StringIO
import itertools


class Hit(object):
    def __init__(self, keyword, rank, database_record_id, fields):
        self.keyword = keyword # string
        self.rank = rank # int
        self.fields = fields # [string]
        self.database_record_id = database_record_id # int

    def __str__(self):
        return '(kw: {}, rank: {}, ncbi id: {}, fields: {})'.format(self.keyword, self.rank, self.database_record_id, ';'.join(self.fields))


class GeneParser(object):
    
    @staticmethod
    def keyword_hit(xmlfile, keywords):
        # parse given xmlfile and extract what we are interested in.
        # look up keywords from extracted content of gene.
        parser = etree.iterparse(xmlfile, events=('end',), tag='Entrezgene')
        for event, geneEntry in parser:    
            geneId, content = get_geneContent(geneEntry)
            for hit in search_keywords_inDict(content, keywords, geneId):
                yield hit


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


# assume input is a set
def fetch_records_from_ids(ids):

    db = 'gene'

    # Before posting, make a list of unique Ids
    #uniqIds = make_unique_list(ids)
    idString = ','.join(list(ids))
    postRequest = Entrez.epost(db=db, id=idString)
    postResult = Entrez.read(postRequest)
    webEnv = postResult['WebEnv']
    queryKey = postResult['QueryKey']
    fetchRequest = Entrez.efetch(db=db,
                                 webenv=webEnv,
                                 query_key=queryKey,
                                 retmode='xml')
    return fetchRequest.read()

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
    geneId = None

    for elem in geneEntry.iterchildren():

        if elem.tag == 'Entrezgene_track-info':
            geneId = elem.find('.//Gene-track/Gene-track_geneid').text
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
    return (geneId, geneContent)



#    def get_keyword_hits_from_pubmed(pmIds):
#        '''Return pubmed_hits'''
#        pubmed_hit = {}
#        for pubmed_content in lookup_pubmed_ids(pmIds):
#            keywords_found = search_keywords_inPubMed(pubmed_content, keywords)
#            pubmed_hit[pubmed_content['PMID']] = keywords_found
#        return pubmed_hit
#
#    def get_keyword_hits_from_genexml(genefilename):
#        '''Return gene_searched and total PMIDs'''
#        genes_searched = {}
#        pmIds = set()
#        geneNames = []
#        # Read genes and extract gene names which we want to find.
#        for name, func, info in read_gene_file(genefilename):
#            if func != 'intergenic':
#                geneNames.append(name)
#
#        content = etree.iterparse(genexmlfile, events=('end',), tag='Entrezgene')
#        for event, geneEntry in content:
#            # Get gene name of geneEntry, and
#            # if the gene name is on the list of gene names which
#            # are in gene file, parse the gene content from xml and
#            # search keywords over the parsed content.
#            geneName = get_geneName(geneEntry)
#            if geneName in geneNames:
#                dbEntry = get_geneContent(geneEntry)
#                searchResult = search_keywords_inDict(dbEntry, keywords)
#                if geneName in genes_searched:
#                    record = genes_searched[geneName]
#                    record['keywordHit'] += searchResult
#                    # XXX this should be a set
#                    record['PmIds'] += make_unique_list(dbEntry['PmIds'])
#                else:
#                    results = {}
#                    results['keywordHit'] = searchResult
#                    results['PmIds'] = make_unique_list(dbEntry['PmIds'])
#                    genes_searched[geneName] = results
#                pmIds.update(results['PmIds'])
#            # Clear node references
#            geneEntry.clear()
#            while geneEntry.getprevious() is not None:
#                del geneEntry.getparent()[0]
#        del content
#        return genes_searched, pmIds
#
#
#    def merge_gene_and_pubmed_hits(genes_searched, pubmed_hit):
#        '''Return genes_searched including pubmed_hits'''
#        # Match gene and pubmed
#        for gene, info_searched in genes_searched.iteritems():
#            gene_pubmed_hit = {}
#            gene_keyword_hit = info_searched['keywordHit']
#            gene_pmids = info_searched['PmIds']
#            # Find pubmed hit information over all PMIDs related to the gene.
#            # Count the number of hits of the keyword over all PMIDs.
#            for pmid in gene_pmids:
#                try:
#                    for keyword, rank in pubmed_hit[pmid]:
#                        if keyword in gene_pubmed_hit:
#                            rank, count = gene_pubmed_hit[keyword]
#                            gene_pubmed_hit[keyword] = (rank, count+1)
#                        else:
#                            gene_pubmed_hit[keyword] = (rank, 1)
#                except KeyError:
#                    pass
#            # Merge keyword sections if the keyword was also found in PMID.
#            for n, (keyword, rank, fields) in enumerate(gene_keyword_hit):
#                if keyword in gene_pubmed_hit:
#                    rank, count = gene_pubmed_hit[keyword]
#                    fields.append('PMID(%s/%s)' % (count, len(gene_pmids)))
#                    gene_keyword_hit[n] = (keyword, rank, fields)
#                    del gene_pubmed_hit[keyword]
#            # For the remaining PMID hits, append to keywordHit.
#            for keyword, (rank, count) in gene_pubmed_hit.items():
#                fields = ['PMID(%s/%s)' % (count, len(gene_pmids))]
#                gene_keyword_hit.append((keyword, rank, fields))
#            # Sort keywordHit according to its rank,
#            # so we can determine top n rank.
#            gene_keyword_hit.sort(key=lambda tup: tup[1])
#            genes_searched[gene]['keywordHit'] = gene_keyword_hit
#        return genes_searched
#

    #genes_searched, pmIds = get_keyword_hits_from_genexml(genefilename)
    #pubmed_hit = get_keyword_hits_from_pubmed(pmIds)
    #genes_searched = merge_gene_and_pubmed_hits(genes_searched, pubmed_hit)
    #return genes_searched




#def search_genes(genefilename, keywords, geneDict):
#    '''Search gene information from geneDict, and
#       print the results.
#           genes: csv.DictReader, keywords: list, geneDict: dict
#       1. searches keywords, and appends results to the original data
#          ['Keyword,Rank']
#    '''
#    column_headers = get_column_headers(genefilename)
#    if column_headers:
#        column_headers.append('Keyword,Rank')
#        print '\t'.join(column_headers)
#
#    for name, func, geneInfo in read_gene_file(genefilename, padding=True):
#        if func == 'intergenic':
#            geneInfo.append('')
#            continue
#        # Search keywords against gene database.
#        result = ''
#        try:
#            dbEntry = geneDict[name]
#            searchResults = search_keywords_inString(dbEntry, keywords)
#            if searchResults:
#                for key, rank in searchResults:
#                    result += '(%s,%s)' % (key, rank)
#            else:
#                # XXX need detail information?
#                result = ''
#        except KeyError:
#            # If DB does not contain the information on gene name.
#            result = "Could not find '%s' in DB" % name
#            print(geneDict.keys())
#            exit()
#        geneInfo.append(result)
#        print '\t'.join(geneInfo)


def search_keywords_inDict(dbEntry, keywords, geneId):
    '''Search top N ranking keywords from dbEntry.
       dbEntry is a dictionary which contains information of each field.
       The value of each key is a list.
       e.g) {'AlterName' : ['abc', 'def'], ...}
       Returns a list of tuples containging
                         the keyword hitted,
                         the rank of keyword, and
                         the fields where the keyword hitted.
    '''
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
            yield Hit(keyword, n+1, geneId, fields)
            fields = []


def search_keywords_inPubMed(pubmedContent, keywords):
    '''Search keywords from pubmedContent.'''
    # Only interested in article title and abstract
    content = ''
    try:
        content += pubmedContent['Article Title']
        content += ';' + pubmedContent['Abstract']
    except KeyError:
        pass

    return search_keywords_inString(content, keywords)


def search_keywords_inString(dbEntry, keywords):
    '''Searche keyword in the order (according to the rank) and
       return the keyword and its rank (position in the list).
       If fail to search, (python) returns None.
       dbEntry is a string.
    '''
    keysFound = []
    for n, item in enumerate(keywords):
        if item in dbEntry:
            keysFound.append((item, n+1))
    return keysFound


def parse_args():
    parser = ArgumentParser(description='Search NCBI for genes of interest, '
                                        'based on concept-keyword search.')

    parser.add_argument('--online',
                        action='store_true',
                        help='Search gene information from online (NCBI).')

    parser.add_argument('--genecol',
                        metavar='INT',
                        type=int,
                        default=0,
                        help='The position of the column containing gene name.')

    parser.add_argument('--skipheader',
                        action='store_true',
                        help='The first line of the gene file is a header which'
                             'should be skipped over')

    parser.add_argument('--organism',
                        type=str,
                        default='human',
                        help='Name of the organism to search')

    parser.add_argument('--email',
                        metavar='EMAIL_ADDRESS',
                        type=str,
                        help='Your email address. This is required by'
                             'NCBI for online queries. You do not need'
                             'to supply an email address for cached queries.')

    parser.add_argument('--genecache',
                        metavar='DIR',
                        type=str,
                        default='genecache',
                        help='Save a cache of the downloaded results '
                             'from NCBI gene into this directory')

    parser.add_argument('--pubmedcache',
                        metavar='DIR',
                        type=str,
                        default='pubmedcache',
                        help='Save a cache of the downloaded results '
                             'from NCBI pubmed into this directory')

    parser.add_argument('--keys',
                        metavar='FILE',
                        type=str,
                        required=True,
                        help='The tab separated file containing '
                             'the keywords to be searched.')

    parser.add_argument('--genes',
                        metavar='FILE',
                        type=str,
                        required=True,
                        help='The tab separated file containing '
                             'the gene information including '
                             'name of the gene, one gene name per line.')

    return parser.parse_args()

