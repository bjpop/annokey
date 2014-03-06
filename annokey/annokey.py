#!/bin/env python

'''
Annokey: a NCBI Gene Database Keyword Search Tool
------------------------------------------------

Authors:   Daniel Park, Sori Kang, Bernie Pope, Tu Nguyen-Dumont.
Copyright: 2013, 2014
Website:   https://github.com/bjpop/annokey
License:   BSD, see LICENCE file in source distribution.

Searches the NCBI gene database for a given set of genes
and a given set of search terms.
'''

import sys
import os
from argparse import ArgumentParser
from Bio import Entrez
import logging
from .version import ANNOKEY_VERSION
#from .report import (DEFAULT_REPORT_FILE, init_report_page, write_report)
from .report import (DEFAULT_REPORT_FILE, report_head, report_foot,
    report_meta_info)
import socket
from .name import PROGRAM_NAME
from .search import search_terms
from .genecache import save_gene_cache

DEFAULT_LOG_FILE = 'annokey_log.txt'

# character used to separate fields in the genes input file.
# Can be a comma or a tab.
DEFAULT_GENE_DELIMITER = ','


def parse_args():
    'A parser for the command line arguments'

    parser = ArgumentParser(description='Search NCBI for genes of interest, '
                                        'based on concept-keyword search.')

    parser.add_argument(
        '--version', action='version', version='%(prog)s ' + ANNOKEY_VERSION)

    parser.add_argument('--online',
                        action='store_true',
                        help='Search gene information from online (NCBI).')

    parser.add_argument('--cachesnapshot',
                        metavar='FILE',
                        type=str,
                        help='Populate the gene cache from downloaded XML '
                              'snapshot of NCBI gene database')

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
                        help='log progress in FILENAME, defaults to {}'.
                             format(DEFAULT_LOG_FILE),
                        default=DEFAULT_LOG_FILE)

    parser.add_argument('--delimiter', type=str,
                        choices=['comma', 'tab'],
                        help='Delimiter for gene file.')

    parser.add_argument('--report', metavar='FILENAME', type=str,
                        help='Save a detailed search report as HTML page, '
                             'defaults to {}'.format(DEFAULT_REPORT_FILE),
                        default=DEFAULT_REPORT_FILE)

    parser.add_argument('--allmatches',
                        help='Return all the matches of a search term in a '
                             'database field, not just the first one',
                        action='store_true')

    parser.add_argument('--pubmed',
                       help='Search titles and abstracts in Pubmed entries '
                            'referred to in each gene entry',
                       action='store_true')

    return parser


def get_gene_delimiter(args):
    '''Determine the kind of delimiter to use in the Genes file.

    1. Check if the --delimiter flag was set to something
    sensible.
    2. Otherwise, check if the filename ends in .csv or .tsv
    3. Otherwise, use the default delimiter and emit a log
       message to that effect.
    '''

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
    '''Entry point for the program.'''

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
    logging.info('annokey version: {}'.format(ANNOKEY_VERSION))
    logging.info('hostname: {}'.format(hostname))
    logging.info('working directory: {}'.format(working_directory))
    logging.info('command line: {}'.format(command_line_text))

    if not ((args.terms and args.genes) or args.cachesnapshot):
        print('\nERROR: Annokey requires --terms AND --genes OR '
              '--cachesnapshot\n')
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
                exit('{}: an email address is required for online/pubmed '
                     'queries, use the --email flag'.format(PROGRAM_NAME))

        args.delimiter = get_gene_delimiter(args)

        with open(args.report, "w") as report_file:
            report_head(report_file)
            report_meta_info(report_file, ANNOKEY_VERSION, hostname,
                working_directory, command_line_text)
            search_terms(args, report_file)
            report_foot(report_file)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        exit('{}: interrupted'.format(PROGRAM_NAME))
