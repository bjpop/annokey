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

import os
from Bio import Entrez
from lxml import etree
from StringIO import StringIO
import itertools


class Hit(object):
    def __init__(self, search_term, rank, database_record_id, fields):
        self.search_term = search_term # string
        self.rank = rank # int
        self.fields = fields # [string]
        self.database_record_id = database_record_id # int

    def __str__(self):
        return '(kw: {}, rank: {}, ncbi id: {}, fields: {})'.format(
            self.search_term,self.rank, self.database_record_id, ';'.join(self.fields))


class GeneParser(object):
    
    @staticmethod
    def term_hit(xmlfile, search_terms, pubmed_cachedir):
        # parse given xmlfile and extract what we are interested in.
        # look up search_terms from extracted content of gene.
        parser = etree.iterparse(xmlfile, events=('end',), tag='Entrezgene')
        for event, geneEntry in parser:
            # XXX We use intermidiate Dictionary at the moment rather than
            # scanning XML file directly, since iterating dictionary N times
            # may be faster than scanning XML file N times when
            # there are N search_terms.
            geneId, content = get_geneContent(geneEntry)
            pubmed_hits = search_terms_in_pubmed(pubmed_cachedir, content["PmIds"], search_terms, geneId)
            for hit in search_terms_inDict(content, search_terms, geneId):
                # If search_term is also in PubMed, append PubMed to field.
                if hit.search_term in pubmed_hits:
                    hit.fields += pubmed_hits[hit.search_term].fields
                    del pubmed_hits[hit.search_term]
                yield hit
            # For the remaining PubMed hits
            for search_term, hit in pubmed_hits.iteritems():
                yield hit

    @staticmethod
    def pubmed_ids(xmlfile):
        # parse given xmlfile and extract pubmed ids.
        pmids = []
        parser = etree.iterparse(xmlfile, events=('end',), tag='Entrezgene')
        for event, geneEntry in parser:
            elem = geneEntry.find('.//Entrezgene_comments/Gene-commentary/Gene-commentary_refs')
            if elem is not None:
                for pub in elem.iterchildren():
                    pubMedId = pub.find('.//Pub_pmid/PubMedId')
                    if pubMedId is not None and pubMedId.text is not None:
                        pmids.append(pubMedId.text)
        return set(pmids)


class PubMedParser(object):

    @staticmethod
    def term_hit(xmlfile, search_terms):
        # Parse the given xmlfile and return search_terms hit.
        parser = etree.iterparse(xmlfile, events=('end',), tag='PubmedArticle')
        for event, pubmed_entry in parser:
            content =  get_pubmed_content(pubmed_entry)
            # scan title and abstract only for search_term search.
            title_abst = ''
            try:
                title_abst += content['Article Title']
                title_abst += content['Abstract']
            except KeyError:
                pass
            termsHit = []
            for search_term in search_terms:
                if search_term.search(title_abst):
                    termsHit.append(search_term)
            return termsHit


# For caching PubMed search_term search history
pubmed_hit_cache = {}

# Search search_terms in PubMed articles. Return a list of Hit.
# We manange pubmed_hit_cache to avoid repeated cache lookup. 
def search_terms_in_pubmed(cachedir, ids, search_terms, gene_id):
    search_term_hits = {}
    for id in ids:
        # For each PubMed article, look up hit history first,
        # If no history in hit cache, look up PubMed file.
        search_term_hit = []
        try:
            search_term_hit = pubmed_hit_cache[id]
        except KeyError:
            pubmed_file = lookup_pubmed_cache(cachedir, id)
            if pubmed_file is not None:
                search_term_hit = PubMedParser.term_hit(pubmed_file, search_terms)
                pubmed_hit_cache[id] = search_term_hit
            else:
                pass
                # XXX maybe should log an error
                #print("did not find PubMed {} in cache:".format(id))

        # Increase hit counts.
        for search_term in search_term_hit:
            try:
                search_term_hits[search_term] += 1
            except KeyError:
                search_term_hits[search_term] = 1

    # Make a list of Hit
    for n, search_term in enumerate(search_terms):
        # XXX there are duplicated search_terms with different ranks.
        # So, we need to check the count is int or not.
        if search_term in search_term_hits:
            count = search_term_hits[search_term]
            if isinstance(count, int):
                field = ["PMID({}/{})".format(search_term_hits[search_term], len(ids))]
                search_term_hits[search_term] = Hit(search_term, n+1, gene_id, field)
    return search_term_hits


def lookup_pubmed_cache(cachedir, id):
    hashed_id = hash(id) % 256
    pubmed_filename = os.path.join(cachedir, str(hashed_id), id)
    if os.path.isfile(pubmed_filename):
        return pubmed_filename


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
            #print("did not find {} in cache:".format(id))
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

    for term, value in values:
        geneContent[term] += value 
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

# XXX
# Dictionary vs Class having attributes
# Not sure about the benefits of Class compared to Dictionary. 
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


def search_terms_inDict(dbEntry, search_terms, geneId):
    '''Search top N ranking search_terms from dbEntry.
       dbEntry is a dictionary which contains information of each field.
       The value of each key is a list.
       e.g) {'AlterName' : ['abc', 'def'], ...}
       Returns a list of tuples containging
                         the search_term hitted,
                         the rank of search_term, and
                         the fields where the search_term hitted.
    '''
    fields = []
    seen_fields = set()
    for n, search_term in enumerate(search_terms):
        for field, content in dbEntry.iteritems():
            if field == 'GeneRIFs':
                hit = [search_term.search(item) for item in content]
                if sum(hit) > 0:
                    fields.append('%s(%s/%s)' %
                                  (field, sum(hit), len(content)))
            else:
                for item in content:
                    if search_term.search(item) and field not in fields:
                        fields.append(field)
        if len(fields) > 0:
            yield Hit(search_term, n+1, geneId, fields)
            fields = []


def search_terms_inString(dbEntry, search_terms):
    '''Searche search_term in the order (according to the rank) and
       return the search_term and its rank (position in the list).
       If fail to search, (python) returns None.
       dbEntry is a string.
    '''
    termsFound = []
    for n, item in enumerate(search_terms):
        if item in dbEntry:
            termsFound.append((item, n+1))
    return termsFound
