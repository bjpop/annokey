'''
Process the XML structure of the NCBI gene database and
the Pubmed database.
'''

import os
from Bio import Entrez
from lxml import etree
from StringIO import StringIO
import logging
from .pubmedcache import (make_pubmed_cache_dirname, save_pubmed_cache)

# XXX maybe should be a named tuple
class GeneHit(object):
    '''A hit for a key term in a field.
       Records:
           - the search term that was used
           - the rank of the search term
           - the name of the database field where the match was made
           - the match object
    '''
    def __init__(self, search_term, rank, field, match):
        self.search_term = search_term # string
        self.rank = rank # int
        self.field = field # string, name of location where hit occurred
        self.match = match # TermMatch

class GeneParser(object):
    '''Processing a Gene XML record, and associated pubmed ids'''
    @staticmethod
    def term_hit(args, xml, search_terms):
        '''For a given gene XML entry find all the hits for all
        the given search terms. Yield the hits one at a time.
        '''
        # parse given xml and extract what we are interested in.
        # look up search_terms from extracted content of gene.
        parser = etree.iterparse(StringIO(xml), events=('end',),
            tag='Entrezgene')
        for _, gene_entry in parser:
            gene_id, content = get_gene_content(gene_entry)
            for hit in search_terms_in_dict(content, search_terms):
                yield hit
            if args.pubmed:
                for hit in search_terms_in_pubmed(args, content["PmIds"],
                    search_terms):
                    yield hit

    @staticmethod
    def pubmed_ids(xmlfile):
        '''Parse given Gene xmlfile and extract pubmed ids'''
        pmids = []
        parser = etree.iterparse(xmlfile, events=('end',), tag='Entrezgene')
        for _, gene_entry in parser:
            elem = gene_entry.find('.//Entrezgene_comments/'
                                  'Gene-commentary/Gene-commentary_refs')
            if elem is not None:
                for pub in elem.iterchildren():
                    pubmed_id = pub.find('.//Pub_pmid/PubMedId')
                    if pubmed_id is not None and pubmed_id.text is not None:
                        pmids.append(pubmed_id.text)
        return set(pmids)


class PubMedParser(object):
    '''Parse pubmed article and find hits for search terms'''

    @staticmethod
    def term_hit(xml, search_terms):
        '''For the XML content of a pubmed article, search for
        each of the search terms. Yield each hit one at a time.
        '''
        # Parse the given xmlfile and return search_terms hit.
        parser = etree.iterparse(StringIO(xml), events=('end',),
            tag='PubmedArticle')
        for _, pubmed_entry in parser:
            content = get_pubmed_content(pubmed_entry)
            # scan title and abstract only for search_term search.
            title_abst = ''
            try:
                title_abst = content['Article Title'] + ' ' + \
                    content['Abstract']
            except KeyError:
                pass
            for rank, search_term in enumerate(search_terms, 1):
                match = search_term.search(title_abst)
                if match is not None:
                    yield GeneHit(search_term, rank, 'PubMed', match)

def search_terms_in_pubmed(args, ids, search_terms):
    '''Lookup pubmed records by their IDs and then search for each
    term in the given list of terms. Yield each hit one at a time.
    '''
    for pubmed_record in get_pubmed_records(args, ids):
        for hit in PubMedParser.term_hit(pubmed_record, search_terms):
            yield hit

# Keep a in-memory cache of pubmed records indexed by the Pubmed id.
# Stops us from requesting it again from online, or from looking again
# in the file cache.
SEEN_PUBMED_RECORDS = {}

def get_pubmed_records(args, pubmed_ids):
    '''Retrieve the pubmed XML records for the given list of pubmed ids.
    We look for the records in the file cache first, before trying to
    fetch them online.
    '''

    not_seen_ids = set()

    # Check for records that we've previously already retrieved.
    for pubmed_id in pubmed_ids:
        if pubmed_id in SEEN_PUBMED_RECORDS:
            yield SEEN_PUBMED_RECORDS[pubmed_id]
        else:
            not_seen_ids.add(pubmed_id)

    not_filecached_ids = set()

    # look for the records in the file cache
    for pubmed_id in not_seen_ids:
        pubmed_filename = lookup_pubmed_cache(args.pubmedcache, pubmed_id)
        if pubmed_filename is not None:
            try:
                with open(pubmed_filename) as pubmed_file:
                    pubmed_xml = pubmed_file.read()
                    SEEN_PUBMED_RECORDS[pubmed_id] = pubmed_xml
                    yield pubmed_xml
            except EnvironmentError:
                logging.warn("Could not open or read pubmed file: {}"
                    .format(pubmed_filename))
        else:
            logging.info("Could not find pubmed article {} in cache"
                .format(pubmed_id))
            not_filecached_ids.add(pubmed_id)

    # Each record from online search could contain many pubmed records.
    # So we split them up and save them individually in the
    # SEEN_PUBMED_RECORDS dictionary.
    for record in fetch_pubmed_records_online(not_filecached_ids):
        parser = etree.iterparse(StringIO(record), events=('end',),
            tag='PubmedArticle')
        for _, pubmed_entry in parser:
            # Find pubmed id
            pubmed_id = pubmed_entry.find('.//MedlineCitation/PMID').text
            logging.info("Found pubmed article {} online".format(pubmed_id))
            record = etree.tostring(pubmed_entry)
            SEEN_PUBMED_RECORDS[pubmed_id] = record
            save_pubmed_cache(args.pubmedcache, pubmed_id, record)
            # free up memory used by the XML iterative parser
            pubmed_entry.clear()
            while pubmed_entry.getprevious() is not None:
                del pubmed_entry.getparent()[0]
            yield record


def fetch_pubmed_records_online(ids):
    '''Fetch the pubmed records from the online Pubmed database
    given a list of their Pubmed Ids.
    '''
    if len(ids) > 0:
        id_string = ','.join(ids)
        post_request = Entrez.epost(db='pubmed', id=id_string)
        post_result = Entrez.read(post_request)
        web_env = post_result['WebEnv']
        query_key = post_result['QueryKey']
        retstart = 0
        chunk_size = 10000

        while retstart < len(ids):
            try:
                fetch_request = Entrez.efetch(db='pubmed',
                                              webenv=web_env,
                                              query_key=query_key,
                                              retmode='xml',
                                              retmax=chunk_size,
                                              retstart=retstart)
                yield fetch_request.read()

            except Exception as exception:
                logging.warn("pubmed fetch failed: {}".format(exception))
                break

            retstart += chunk_size


def lookup_pubmed_cache(cachedir, pubmed_id):
    '''Try to find the file containing a given Pubmed article
    from the Pubmed file cache based on its Pubmed id.
    If it can't be found return None
    '''
    pubmed_cache_dir = make_pubmed_cache_dirname(cachedir, pubmed_id)
    pubmed_cache_filename = os.path.join(pubmed_cache_dir, pubmed_id)
    if os.path.isfile(pubmed_cache_filename):
        return pubmed_cache_filename


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


def get_pubmed_content(pubmed_entry):
    '''Parse element PubmedArticle'''
    pubmed_content = {}
    for elem in pubmed_entry.iterchildren():
        if elem.tag == 'MedlineCitation':
            pmid = elem.find('.//PMID')
            if pmid is not None and pmid.text is not None:
                pubmed_content['PMID'] = pmid.text
            article = elem.find('.//Article')
            if article is not None:
                journal = elem.find('.//Journal')
                article_title = elem.find('.//ArticleTitle')
                pagination = elem.find('.//Pagination')
                abstract = elem.find('.//Abstract')
            if journal is not None:
                volume = journal.find('.//JournalIssue/Volume')
                issue = journal.find('.//JournalIssue/Issue')
                pubdate = journal.find('.//JournalIssue/PubDate/Year')
                journal_title = journal.find('.//Title')
            if volume is not None and volume.text is not None:
                pubmed_content['Volume'] = volume.text
            if issue is not None and issue.text is not None:
                pubmed_content['Issue'] = issue.text
            if pubdate is not None and pubdate.text is not None:
                pubmed_content['Year'] = pubdate.text
            if journal_title is not None and journal_title.text is not None:
                pubmed_content['Journal Title'] = journal_title.text
            if article_title is not None and article_title.text is not None:
                pubmed_content['Article Title'] = article_title.text
            if pagination is not None and pagination.text is not None:
                pgn = pagination.find('.//MedlinePgn')
                if pgn is not None and pgn.text is not None:
                    pubmed_content['Pagination'] = pgn.text
            if abstract is not None:
                abstract_text = ''
                for abstract_child in abstract.iterchildren():
                    if abstract_child.tag == 'AbstractText':
                        if abstract_child.text is not None:
                            abstract_text += abstract_child.text
                pubmed_content['Abstract'] = abstract_text
    return pubmed_content


def merge_gene_content(gene_content, values):
    '''Merge gene information.

    Args:
         gene_content: dict containing gene information.
                      values to be merged into this gene_content.
         values: a list of tuple containing gene information.
                 It is to be merged into gene_content.
    Returns:
         A dict containing gene information.
    '''

    for term, value in values:
        gene_content[term] += value
    return gene_content


def get_other_source_anchor(entry):
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


def parse_gene_commentary(commentary):
    '''Parse element 'Gene-commentary'''

    return_values = []

    rifs = []
    pathways = []
    pmids = set()
    interactions = []

    element_type = ''
    heading = ''
    label = ''

    for item in commentary.iterchildren():
        if item.tag == 'Gene-commentary_type':
            element_type = item.attrib['value']

        if item.tag == 'Gene-commentary_heading':
            heading = item.text

        if item.tag == 'Gene-commentary_label':
            label = item.text

        if item.tag == 'Gene-commentary_text':
            if element_type == 'generif' and item.text is not None:
                rifs.append(item.text)
                return_values.append(('GeneRIFs', rifs))

        if item.tag == 'Gene-commentary_refs':
            for pub in item.iterchildren():
                pubmed_id = pub.find('.//Pub_pmid/PubMedId')
                if pubmed_id is not None and pubmed_id.text is not None:
                    pmids.add(pubmed_id.text)
            return_values.append(('PmIds', pmids))

        if item.tag == 'Gene-commentary_products':
            for child in item.iterchildren():
                return_values += parse_gene_commentary(child)

        if item.tag == 'Gene-commentary_comment':
            if heading == 'Pathways':
                for child in item.iterchildren():
                    pathway = child.find('.//Gene-commentary_text')
                    if pathway is not None and pathway.text is not None:
                        pathways.append(pathway.text)
                return_values.append(('Pathways', pathways))

            elif heading == 'Interactions':
                interactions = []
                for child in item.iterchildren():
                    #sub_comment = item.find('.//Gene-commentary_comment')
                    for grandchild in child.iterchildren():
                        interactions += get_other_source_anchor(grandchild)
                return_values.append(('Interactions', interactions))

            elif heading == 'NCBI Reference Sequences (RefSeq)':
                for child in item.iterchildren():
                    return_values += parse_gene_commentary(child)

            elif heading == 'Conserved Domains':
                result = get_other_source_anchor(item)
                return_values.append(('Conserved Domains', result))

            elif label == 'Function':
                result = get_other_source_anchor(item)
                return_values.append(('Function', result))

            elif label == 'Component':
                result = get_other_source_anchor(item)
                return_values.append(('Component', result))

            elif label == 'Process':
                result = get_other_source_anchor(item)
                return_values.append(('Process', result))

            else:
                # Find 'Conserved Domains'
                for child in item.iterchildren():
                    if child.tag == 'Gene-commentary':
                        heading = child.find('.//Gene-commentary_heading')
                        if (heading is not None and
                                heading.text == 'Conserved Domains'):
                            return_values += parse_gene_commentary(child)

    return return_values

def get_gene_content(gene_entry):
    '''Parse element Entrezgene.'''

    gene_content = {'Gene Name': [],
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
    gene_id = None

    for elem in gene_entry.iterchildren():

        if elem.tag == 'Entrezgene_track-info':
            gene_id = elem.find('.//Gene-track/Gene-track_geneid').text
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
                gene_content['Gene Name'].append(name.text)
            if desc is not None and desc.text is not None:
                gene_content['Description'].append(desc.text)
            if syn is not None:
                for child in syn.iterchildren():
                    if (child.tag == 'Gene-ref_syn_E' and
                            child.text is not None):
                        gene_content['Synonyms'].append(child.text)

        elif elem.tag == 'Entrezgene_summary' and elem.text is not None:
            gene_content['Summary'].append(elem.text)

        elif elem.tag == 'Entrezgene_comments':
            for sub_elem in elem.iterchildren():
                result = parse_gene_commentary(sub_elem)
                gene_content = merge_gene_content(gene_content, result)

        elif elem.tag == 'Entrezgene_properties':
            entry = elem.find('.//Gene-commentary/Gene-commentary_comment')
            if entry is not None:
                for child in entry.iterchildren():
                    result = parse_gene_commentary(child)
                    gene_content = merge_gene_content(gene_content, result)

        elif elem.tag == 'Entrezgene_prot':
            alter_name = elem.find('.//Prot-ref/Prot-ref_name')
            if alter_name is not None:
                for name in alter_name.iterchildren():
                    if (name.tag == 'Prot-ref_name_E' and
                            name.text is not None):
                        gene_content['Alternative Name'].append(name.text)
    return (gene_id, gene_content)


def search_terms_in_dict(db_entry, search_terms):
    '''Iteratate over all the entries in a dictionary and search
    for each search term in each entry.
    '''
    for rank, search_term in enumerate(search_terms, 1):
        for field, content in db_entry.iteritems():
            for item in content:
                match = search_term.search(item)
                if match is not None:
                    yield GeneHit(search_term, rank, field, match)
