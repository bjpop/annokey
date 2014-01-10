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
from version import annokey_version
from genecache import (lookup_gene_cache_iter, stable_string_hash,
                      make_gene_cache_dirname, save_gene_cache)

DEFAULT_LOG_FILE = 'annokey_log.txt'

# character used to separate fields in the genes input file.
# Can be a comma or a tab.
DEFAULT_GENE_DELIMITER = ','


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


# read each gene name from the gene file
def get_gene_names(args, genefilename):
    with open(genefilename, 'rb') as file:
        reader = csv.DictReader(file, delimiter=args.delimiter)
        for row in reader:
            name = row['Gene'].strip()
            yield name


def get_gene_ids(args, genefilename, organism='Homo sapiens'):

    def chunk_gene_names(genefilename, chunk_size):
        names = []
        numbered_gene_names = enumerate(get_gene_names(args, genefilename))
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


# For each gene in genesfile we read PubMedIds using genecache.
def get_pubmed_ids(args):
    pubmed_ids = set()
    
    with open(args.genes) as genesfile:
        reader = csv.DictReader(genesfile, delimiter=args.delimiter)
        for row in reader:
            genename = row['Gene'].strip()
            for gene_xml in lookup_gene_cache_iter(args, genename):
                pubmed_ids.update(GeneParser.pubmed_ids(gene_xml))

    return pubmed_ids


def fetch_and_save_pubmed_records(cachedir, ids):
    # XXX Fetch all records here regardless of existing pubmed cache.
    # We may be able to skip the records that already exist if
    # the records are not updated.
    # XXX As the number of pubmed ids is usually large,
    # save XML file in cache right after fetching them.
    if len(ids) > 0:
        idString = ','.join(ids)
        postRequest = Entrez.epost(db='pubmed', id=idString)
        postResult = Entrez.read(postRequest)
        webEnv = postResult['WebEnv']
        queryKey = postResult['QueryKey']
        retstart = 0
        chunk_size = 10000

        while retstart < len(ids):
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

            save_pubmed_cache(cachedir, pubmed_result)
            retstart += chunk_size


def lookup_pubmed_ids(ids):
    not_cached_ids = []
    # search for all the cached pubmed ids first, and
    # collect the non-cached ones.
    for id in ids:
        cache_result = lookup_pubmed_cache(id)
        if cache_result is not None:
            yield cache_result
        else:
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


def summarise_hits(hits):
    '''Given a list of hits from a keyword search, generate a list with
    the following three things:

      1. highest ranked hit
      2. all the terms that were in hits
      3. all the fields where the terms were found
    '''

    ranks = []
    terms = set()
    fields = set()

    for h in hits:
        ranks.append(h.rank)
        terms.add(h.keyword)
        for f in h.fields:
            fields.add(f)

    if len(ranks) > 0:
        highest_rank = str(sorted(ranks)[0])
    else:
        highest_rank = ''

    all_terms = ';'.join([t.strip() for t in terms])
    all_fields = ';'.join([f.strip() for f in fields])

    return [highest_rank, all_terms, all_fields]


def search_keywords(args):

    # build a list of all the keywords in the order that they
    # appear in the keywords file (rank order)
    keywords = []
    with open(args.keys) as keysfile:
        for line in keysfile:
            keywords.append(line.strip())

    if len(keywords) > 0:
        with open(args.genes) as genesfile:
            reader = csv.DictReader(genesfile, delimiter=args.delimiter)
            writer = csv.writer(sys.stdout, delimiter=args.delimiter)
            # preserve the header from the original gene file
            annokey_headers = ['Highest Ranked Match', 'Matched Terms', 'Fields']
            new_header = reader.fieldnames + annokey_headers 
            writer.writerow(new_header)
            # for each row in the genes file, find hits for the gene
            # print the row out with the hits annoated on the end
            for input_row in reader:
                try:
                    genename = input_row['Gene'].strip()
                except KeyError:
                    exit("Can't find Gene column in input gene file")
                else:
                    # hits = [str(hit) for hit in  search_keywords_gene_iter(args, genename, keywords)] 
                    hits = list(search_keywords_gene_iter(args, genename, keywords))
                    hits_output = summarise_hits(hits) 
                    output_row = [input_row[field] for field in reader.fieldnames]
                    writer.writerow(output_row + hits_output)


# Search for each keyword in the XML file for a gene. A hit is yielded for
# each keyword. If a gene has multiple files, we search each one separately.
# This means it is possible to get multiple hits for the same keyword
# but from different files. XXX we should annotate each hit with the database
# file ID
def search_keywords_gene_iter(args, gene_name, keywords):
    for gene_xml in lookup_gene_cache_iter(args, gene_name):
        for hit in GeneParser.keyword_hit(gene_xml, keywords, args.pubmedcache):
            yield hit


def make_pubmed_cache_dirname(cachedir, pubmed_id):
    hash_dir = stable_string_hash(pubmed_id) % 256
    return os.path.join(cachedir, str(hash_dir))


def save_pubmed_cache(cachedir, pubmed_records_xml):
    # Read each "PubMedArticle" record in the input XML and 
    # write it out to a cache file. We store each article entry in a file
    # based on its PubMed ID.
    
    parser = etree.iterparse(StringIO(pubmed_records_xml), events=('end',), tag='PubmedArticle')
    for event, pubmed_entry in parser:
        # Find pubmed id
        pubmed_id = pubmed_entry.find('.//MedlineCitation/PMID').text
        pubmed_cache_dir = make_pubmed_cache_dirname(cachedir, pubmed_id)
        if not os.path.exists(pubmed_cache_dir):
            os.makedirs(pubmed_cache_dir)
        pubmed_cache_filename = os.path.join(pubmed_cache_dir, pubmed_id)
        with open(pubmed_cache_filename, 'w') as cache_file:
            cache_file.write(etree.tostring(pubmed_entry))
        # free up memory used by the XML iterative parser
        pubmed_entry.clear()
        while pubmed_entry.getprevious() is not None:
            del pubmed_entry.getparent()[0]


def fetch_records(args):
    # fetch gene records and pubmed records and save them in cache.

    # read each gene name from the gene file 
    gene_ids = get_gene_ids(args, args.genes, args.organism)
    # fetch the corresponding XML records for the identified genes
    gene_records_xml = fetch_records_from_ids(gene_ids)
    # Cache the gene entries in the XML file
    save_gene_cache(args.genecache, args.organism, StringIO(gene_records_xml))

    # pubmed records
    # get pubmed ids from each gene records
    pubmed_ids = get_pubmed_ids(args)
    # fetch the corresponding XML records for pubmed ids
    # cache the pubmed records
    fetch_and_save_pubmed_records(args.pubmedcache, pubmed_ids)


def parse_args():
    parser = ArgumentParser(description='Search NCBI for genes of interest, '
                                        'based on concept-keyword search.')

    parser.add_argument(
        '--version', action='version', version='%(prog)s ' + annokey_version)

    parser.add_argument('--online',
                        action='store_true',
                        help='Search gene information from online (NCBI).')

    parser.add_argument('--xml',
                        metavar='FILE',
                        type=str,
                        help='Populate gene information from downloaded XML dump of NCBI gene database')

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

    parser.add_argument('--log', metavar='FILENAME', type=str,
                        help='log progress in FILENAME, defaults to {}'.format(DEFAULT_LOG_FILE),
                        default=DEFAULT_LOG_FILE)

    parser.add_argument('--delimiter', type=str,
                        choices=['comma', 'tab'],
                        help='delimiter for gene file')

    return parser.parse_args()


def get_gene_delimiter(args):
    if args.delimiter == 'comma':
        return ','
    elif args.delimiter == 'tab':
        return '\t'
    elif args.genes.endswith('.csv'):
        return ','
    elif args.genes.endswith('.tsv'):
        return '\t'
    else:
        logging.info("Using '{}' for delimiter in genes input file"
                     .format(DEFAULT_GENE_DELIMITER))
        return DEFAULT_GENE_DELIMITER


def main():
    args = parse_args()

    if args.log:
        logging.basicConfig(filename=args.log,
                        level=logging.DEBUG,
                        filemode='w',
                        format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))

    args.delimiter = get_gene_delimiter(args)

    if args.xml:
        # Get gene information from specified XML file.
        # Populate the genecache from the contents of the file. 
        with open(args.xml) as xml_file:
            save_gene_cache(args.genecache, args.organism, xml_file)

    elif args.online:
        # Get gene information online.
        if args.email:
            Entrez.email = args.email
        else:
            exit('An email address is required for online queries, use the --email flag')

        # fetch gene and pubmed records and save them in cache. 
        fetch_records(args)

    # Get gene information from cache.
    search_keywords(args)

if __name__ == '__main__':
    main()
