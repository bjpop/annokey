#!/bin/env python

'''
--------------------------------------------------------------------------------

Annokey: a NCBI Gene Database Keyword Search Tool
------------------------------------------------

Authors:   Daniel Park, Sori Kang, Bernie Pope, Tu Nguyen-Dumont.
Copyright: 2013
Website:   https://github.com/bjpop/annokey
License:   BSD, see LICENCE file in source distribution. 

Searches the NCBI gene database for a given set of genes
and a given set of search terms. 

This program supports both online search and offline search.
If the option --online is on, this program downloads gene database in xml
format from the NCBI server and search terms. The option --saveCache allows
the downloaded information to be saved int the user's directory.

For offline search, there is the option --loadCache, which
accepts an xml file for gene information. 

The search results are appended at the end of column of the
input gene file with the information about the search term hit,
the search term rank, and the section where the search term hit.

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
import logging
from process_xml import GeneParser, GeneHit
from version import annokey_version
from genecache import (lookup_gene_cache_iter,
    make_gene_cache_dirname, save_gene_cache)
from pubmedcache import make_pubmed_cache_dirname
from report import (DEFAULT_REPORT_FILE, init_report_page, report_hits, write_report)
from search_term import (parse_search_term)
import socket
from name import program_name

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


def aggregate_hits(hits):
    '''Aggregate information about all the hits for a search term
    and a given gene.

    Result is a dictionary indexed by search term.

    Each search term points to a dictionary indexed by field name
    (the place in the database where the term was matched).

    Each field name points to a list of contexts (strings surrounding
    the matched term.
    '''
    term_dict = {}
    for h in hits:
        term = str(h.search_term)
        field = h.field
        match = h.match
        rank = h.rank
        if (rank, term) in term_dict:
            term_fields = term_dict[(rank, term)]
            if field in term_fields:
                term_field_matches = term_fields[field]
                term_field_matches.append(match)
            else:
                term_fields[field] = [match]
        else:
            term_fields = {field: [match]}
            term_dict[(rank, term)] = term_fields
    return term_dict


def summarise_hits(all_hits):

    all_ranks = []
    num_matches = 0

    for (rank, term), fields in all_hits.items():
        all_ranks.append(rank)
        for field, matches in fields.items():
            for match in matches:
                num_matches += len(match.spans)
      
    if len(all_ranks) > 0:
        highest_rank = str(sorted(all_ranks)[0])
        return [highest_rank, str(num_matches)]
    else:
        return ['', '']


def search_terms(args, report_page):

    # build a list of all the search terms in the order that they
    # appear in the search terms file (rank order)
    terms = []

    try:
        with open(args.terms) as termsfile:
            for line in termsfile:
                terms.append(parse_search_term(line))
    except EnvironmentError as e:
        exit("{}: failed to open terms file: {}".format(program_name, e))

    if len(terms) > 0:
        with open(args.genes) as genesfile:
            reader = csv.DictReader(genesfile, delimiter=args.delimiter)
            writer = csv.writer(sys.stdout, delimiter=args.delimiter)
            # preserve the header from the original gene file
            annokey_headers = ['Highest Ranked Match', 'Num Matched Entries']
            new_header = reader.fieldnames + annokey_headers 
            writer.writerow(new_header)
            # for each row in the genes file, find hits for the gene
            # print the row out with the hits annoated on the end
            for input_row in reader:
                try:
                    gene_name = input_row['Gene'].strip()
                except KeyError:
                    exit("{}: can't find Gene column in input gene file".format(program_name))
                else:
                    for gene_db_id, gene_xml in get_gene_records(args, gene_name):
                        hits = list(GeneParser.term_hit(args, gene_xml, terms))
                        all_hits = aggregate_hits(hits)
                        summary  = summarise_hits(all_hits)
                        output_row = [input_row[field] for field in reader.fieldnames]
                        writer.writerow(output_row + summary)
                        report_hits(gene_name, gene_db_id, all_hits, report_page)


def get_gene_records(args, gene_name):
    if args.online:
        gene_ids = get_gene_ids(args, gene_name)
        for gene_id in gene_ids:
            gene_records_xml = fetch_records_from_ids(gene_ids)
            yield gene_id, gene_records_xml
    else:
       for gene_db_id, gene_file_path in lookup_gene_cache_iter(args, gene_name):
           try:
               with open(gene_file_path) as gene_xml_file:
                   gene_xml_contents = gene_xml_file.read()
           except EnvironmentError as e:
               exit("{}: failed to open gene XML file: {}".format(gene_file_path, e))
           yield gene_db_id, gene_xml_contents

def get_gene_ids(args, gene_name):
    genefilename = args.genes
    organism = args.organism
    search_term = '{}[organism] AND {}[sym]'.format(organism, gene_name)
    request = Entrez.esearch(db='gene', term=search_term, retmax=1000)
    try:
        result = Entrez.read(request)['IdList']
    except Exception as e:
        logging.warn("Cannot get gene ids for gene {}: {}".format(gene_name, e))
        return []
    else:
        return result

def fetch_records_from_ids(ids):
    db = 'gene'
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


def parse_args():
    parser = ArgumentParser(description='Search NCBI for genes of interest, '
                                        'based on concept-keyword search.')

    parser.add_argument(
        '--version', action='version', version='%(prog)s ' + annokey_version)

    parser.add_argument('--online',
                        action='store_true',
                        help='Search gene information from online (NCBI).')

    parser.add_argument('--cachesnapshot',
                        metavar='FILE',
                        type=str,
                        help='Populate the gene cache from downloaded XML snapshot of NCBI gene database')

    parser.add_argument('--organism',
                        type=str,
                        default='human',
                        help='Name of the organism to search')

    parser.add_argument('--email',
                        metavar='EMAIL_ADDRESS',
                        type=str,
                        help='Your email address. This is required by '
                             'NCBI for online queries. You do not need '
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

    parser.add_argument('--terms',
                        metavar='FILE',
                        type=str,
                        help='The tab separated file containing '
                             'the search terms to be searched.')

    parser.add_argument('--genes',
                        metavar='FILE',
                        type=str,
                        help='The tab separated file containing '
                             'the gene information including '
                             'name of the gene, one gene name per line.')

    parser.add_argument('--log', metavar='FILENAME', type=str,
                        help='log progress in FILENAME, defaults to {}'.format(DEFAULT_LOG_FILE),
                        default=DEFAULT_LOG_FILE)

    parser.add_argument('--delimiter', type=str,
                        choices=['comma', 'tab'],
                        help='Delimiter for gene file.')

    parser.add_argument('--report', metavar='FILENAME', type=str,
                        help='Save a detailed search report as HTML page, defaults to {}'.format(DEFAULT_REPORT_FILE),
                        default=DEFAULT_REPORT_FILE)

    parser.add_argument('--allmatches',
                        help='Return all the matches of a search term in a database field, not just the first one',
                        action='store_true')

    parser.add_argument('--pubmed',
                       help='Search titles and abstracts in Pubmed entries referred to in each gene entry',
                       action='store_true')

    return parser


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
    args_parser = parse_args()
    args = args_parser.parse_args()

    if args.log:
        logging.basicConfig(filename=args.log,
                        level=logging.DEBUG,
                        filemode='w',
                        format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')
    command_line_text = "annokey " + ' '.join(sys.argv[1:])
    working_directory = os.getcwd()
    hostname = socket.gethostname()
    logging.info('hostname: {}'.format(hostname))
    logging.info('working directory: {}'.format(working_directory))
    logging.info('command line: {}'.format(command_line_text))

    if not ((args.terms and args.genes) or args.cachesnapshot):
        print("\nERROR: Annokey requires --terms AND --genes OR --cachesnapshot\n")
        args_parser.print_help()
        exit()

    if args.cachesnapshot:
        # Get gene information from specified XML file.
        # Populate the genecache from the contents of the file. 
        save_gene_cache(args)

    # only perform search if terms and genes are specified
    if args.genes and args.terms:
        if args.online or args.pubmed:
            #  Get gene information online.
            if args.email:
                Entrez.email = args.email
            else:
                exit('{}: an email address is required for online/pubmed queries, use the --email flag'.format(program_name))

        args.delimiter = get_gene_delimiter(args)
        report_page = init_report_page(hostname, working_directory, command_line_text)
        search_terms(args, report_page)
        write_report(args.report, report_page)


if __name__ == '__main__':
    main()
