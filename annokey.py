#!/bin/env python

'''
--------------------------------------------------------------------------------

Annokey: a NCBI Gene Database Keyword Search Tool
-------------------------------------------------

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

For offline search, there is the option --loadCache, which
accepts an xml file for gene information. 

The search results are appended at the end of column of the
input gene file with the information about the keyword hit,
the keyword rank, and the section where the keyword hit.

Required inputs:

    --keys FILENAME   a text file of search keywords/keyphrases. One
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
import logging

from process_xml import GeneParser, Hit


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


def genefile_reader(csvfile):
    # try to detect the csv format from the (up to)
    # first three lines of the file
    sample = ''
    for n, line in enumerate(csvfile):
        if n >= 3:
            break
        else:
            sample += line 
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(sample)
    csvfile.seek(0)
    reader = csv.reader(csvfile, dialect)
    if sniffer.has_header:
        header = next(reader)
    else:
        header = None
    return header, reader, dialect


# read each gene name from the gene file at a given column
def get_gene_names(genefilename, column):
    with open(genefilename, 'rb') as file:
        header, reader, dialect = genefile_reader(file)
        for row in reader:
            name = row[column]
            yield name


def get_gene_ids(genefilename, column, organism='Homo sapiens'):

    def chunk_gene_names(genefilename, chunk_size):
        names = []
        numbered_gene_names = enumerate(get_gene_names(genefilename, column))
        for n, name in numbered_gene_names:
            names.append(name.strip())
            if (n+1) % chunk_size == 0:
                yield names
                names = []
        yield names

    # Send request to server for every 100 genes.
    gene_ids = set()

    for names in chunk_gene_names(genefilename, chunk_size=100):
        terms = ' OR '.join(['%s[sym]' % name for name in names])
        search_term = '{0}[organism] AND ({1})'.format(organism, terms)
        request = Entrez.esearch(db='gene', term=search_term, retmax=10000)
        result = Entrez.read(request)
        gene_ids.update(result['IdList'])
    return gene_ids


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


#class Hit(object):
#    def __init__(self, keyword, rank, database_record_id, fields):
#        self.keyword = keyword # string
#        self.rank = rank # int
#        self.fields = fields # [string]
#        self.database_record_id = database_record_id # int
#
#    def __str__(self):
#        return '(kw: {}, rank: {}, ncbi id: {}, fields: {})'.format(self.keyword, self.rank, self.database_record_id, ';'.join(self.fields))


def search_keywords(args):

    # build a list of all the keywords in the order that they
    # appear in the keywords file
    keywords = []
    with open(args.keys) as keysfile:
        for line in keysfile:
            keywords.append(line.strip())

    if len(keywords) > 0:
        with open(args.genes) as genesfile:
            header, reader, dialect = genefile_reader(genesfile)
            writer = csv.writer(sys.stdout, dialect=dialect)
            if header:
                writer.writerow(header)
            # for each row in the genes file, find hits for the gene
            # print the row out with the hits annoated on the end
            for row in reader:
                genename = row[args.genecol].strip()
                for hit in search_keywords_gene_iter(args, genename, keywords):
                    row.append(str(hit))
                writer.writerow(row)


# Search for each keyword in the XML file for a gene. A hit is yielded for
# each keyword. If a gene has multiple files, we search each one separately.
# This means it is possible to get multiple hits for the same keyword
# but from different files. XXX we should annotate each hit with the database
# file ID
def search_keywords_gene_iter(args, gene_name, keywords):
    for gene_xml in lookup_gene_cache_iter(args, gene_name):
        for hit in GeneParser.keyword_hit(gene_xml, keywords):
            yield hit


def lookup_gene_cache_iter(args, gene_name):
    gene_cache_dir = make_gene_cache_dirname(args.genecache, args.organism, gene_name)
    try:
        dir_contents = os.listdir(gene_cache_dir)
    except OSError:
        logging.info("Could not find cache entry for {}".format(gene_name))
        return
    for filename in os.listdir(gene_cache_dir):
        #if filename.endswith('.xml'):
        file = '%s/%s' % (gene_cache_dir, filename)
        file = open(file)
        yield file
        file.close()


#def search_keywords_inDict(dbEntry, keywords):
#    '''Search top N ranking keywords from dbEntry.
#       dbEntry is a dictionary which contains information of each field.
#       The value of each key is a list.
#       e.g) {'AlterName' : ['abc', 'def'], ...}
#       Returns a list of tuples containging
#                         the keyword hitted,
#                         the rank of keyword, and
#                         the fields where the keyword hitted.
#    '''
#    keysFound = []
#    fields = []
#    for n, keyword in enumerate(keywords):
#        for field, content in dbEntry.iteritems():
#            if field == 'GeneRIFs':
#                hit = [keyword in item for item in content]
#                if sum(hit) > 0:
#                    fields.append('%s(%s/%s)' %
#                                  (field, sum(hit), len(content)))
#            else:
#                for item in content:
#                    if keyword in item and field not in fields:
#                        fields.append(field)
#        if len(fields) > 0:
#            keysFound.append((keyword, n+1, fields))
#            fields = []
#
#    return keysFound
#
#
#def search_keywords_inPubMed(pubmedContent, keywords):
#    '''Search keywords from pubmedContent.'''
#    # Only interested in article title and abstract
#    content = ''
#    try:
#        content += pubmedContent['Article Title']
#        content += ';' + pubmedContent['Abstract']
#    except KeyError:
#        pass
#
#    return search_keywords_inString(content, keywords)
#
#
#def search_keywords_inString(dbEntry, keywords):
#    '''Searche keyword in the order (according to the rank) and
#       return the keyword and its rank (position in the list).
#       If fail to search, (python) returns None.
#       dbEntry is a string.
#    '''
#    keysFound = []
#    for n, item in enumerate(keywords):
#        if item in dbEntry:
#            keysFound.append((item, n+1))
#    return keysFound



def make_gene_cache_dirname(cachedir, organism, gene_name):
    organism_cache_dir = os.path.join(cachedir, organism, 'gene')
    hash_dir = hash(gene_name) % 256
    return os.path.join(organism_cache_dir, str(hash_dir), gene_name)

def save_gene_cache(cachedir, organism, gene_records_xml):

    #cachedir = args.genecache
    #organism = args.organism

    # create a cache directory for this organism if it doesn't already exist
    # for example human would go in:
    # cachedir/human/gene/
    # the reason for "gene" is because we also cache ncbi entries in "ncbi"
    #organism_cache_dir = os.path.join(cachedir, organism, 'gene')

    # read each "Entrezgene" record in the input XML and write it out to
    # a cache file. Sometimes one gene will have multiple entries. We store
    # each gene entry in a file based on its database ID.
    parser = etree.iterparse(StringIO(gene_records_xml), events=('end',), tag='Entrezgene')

    for event, elem in parser:
        # find the official name of the gene
        gene_name = elem.find('.//Gene-ref_locus').text
        # find the database id of the gene
        gene_id = elem.find('.//Gene-track_geneid').text
        # hash the name of the gene into an integer in the range [0, 255]
        #hash_dir = hash(gene_name) % 256
        # if it doesn't already exist, create a directory in the cache for
        # this gene
        #gene_cache_dir = os.path.join(organism_cache_dir, str(hash_dir), gene_name)
        gene_cache_dir = make_gene_cache_dirname(cachedir, organism, gene_name)
        if not os.path.exists(gene_cache_dir):
            os.makedirs(gene_cache_dir)
        # Write a single 'Entrezgene' entry to a file in the cache using the
        # database ID for the file name
        gene_cache_filename = os.path.join(gene_cache_dir, gene_id)
        with open(gene_cache_filename, 'w') as cache_file:
            cache_file.write(etree.tostring(elem))
        # free up memory used by the XML iterative parser
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]


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
                        required=True,
                        help='The position of the column containing gene name. (0 base)')

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

    parser.add_argument('--log', metavar='FILENAME', type=str, required=True,
                        help='log progress in FILENAME')


    return parser.parse_args()

def main():
    args = parse_args()

    logging.basicConfig(filename=args.log,
                        level=logging.DEBUG,
                        filemode='w',
                        format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))

    if args.online:
        # Get gene information online.

        if args.email:
            Entrez.email = args.email
        else:
            exit('An email address is required for online queries, use the --email flag')

        # read each gene name from the gene file at a given column
        gene_ids = get_gene_ids(args.genes, args.genecol, args.organism)
        # fetch the corresponding XML records for the identified genes
        gene_records_xml = fetch_records_from_ids(gene_ids)

        # Cache the gene entries in the XML file
        save_gene_cache(args.genecache, args.organism, gene_records_xml) 

    # Get gene information from cache.
    search_keywords(args)

if __name__ == '__main__':
    main()
