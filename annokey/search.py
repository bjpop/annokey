'''
Search code for Annokey.

Find matches of key terms within selected fields of
the NCBI gene database and Pubmed articles.
'''

import sys
import csv
from Bio import Entrez
import logging
from .process_xml import GeneParser
from .genecache import lookup_gene_cache_iter
from .report import report_hits
from .name import PROGRAM_NAME
import re

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
    for hit in hits:
        term = str(hit.search_term)
        field = hit.field
        match = hit.match
        rank = hit.rank
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
    '''From the aggregated hits, compute a summary for the purposes of
    annotating the CSV file.

    The annotation contains:

       1) The rank of the highest ranking matched search term.
       2) The highest ranking matched search term (the term itself).
       3) The total number of times the term was matched in all
          the fields.
    '''

    all_ranks = []
    num_matches = 0

    for rank_term, fields in all_hits.items():
        all_ranks.append(rank_term)
        for _, matches in fields.items():
            for match in matches:
                num_matches += len(match.spans)

    if len(all_ranks) > 0:
        highest_rank, highest_rank_term = sorted(all_ranks)[0]
        return [str(highest_rank), highest_rank_term, str(num_matches)]
    else:
        return ['', '', '']


def search_terms(args, report_page):
    '''For each gene in the input, search for occurrences
    of each of the search terms.

    Produce a HTML report of all the hits, and annotate
    and output CSV file with a summary of the hits.
    '''

    # build a list of all the search terms in the order that they
    # appear in the search terms file (rank order)
    terms = []

    try:
        with open(args.terms) as termsfile:
            for line in termsfile:
                terms.append(parse_search_term(line))
    except EnvironmentError as exception:
        exit("{}: failed to open terms file: {}".
            format(PROGRAM_NAME, exception))

    if len(terms) > 0:
        search_terms_worker(args, terms, report_page)

def search_terms_worker(args, terms, report_page):
    '''Worker function for search_terms.'''

    with open(args.genes) as genesfile:
        reader = csv.DictReader(genesfile, delimiter=args.delimiter)
        writer = csv.writer(sys.stdout, delimiter=args.delimiter)
        # preserve the header from the original gene file
        annokey_headers = ['Highest Rank', 'Highest Rank Term',
            'Total Matched Entries']
        new_header = reader.fieldnames + annokey_headers
        writer.writerow(new_header)
        # for each row in the genes file, find hits for the gene
        # print the row out with the hits annoated on the end
        for input_row in reader:
            try:
                gene_name = input_row['Gene'].strip()
            except KeyError:
                exit("{}: can't find Gene column in input gene file".
                    format(PROGRAM_NAME))
            else:
                for gene_count, (gene_db_id, gene_xml) in \
                        enumerate(get_gene_records(args, gene_name)):
                    hits = list(GeneParser.term_hit(args, gene_xml, terms))
                    all_hits = aggregate_hits(hits)
                    summary = summarise_hits(all_hits)
                    output_row = [input_row[field] for field in
                                 reader.fieldnames]
                    writer.writerow(output_row + summary)
                    report_hits(gene_count, gene_name, gene_db_id, all_hits,
                        report_page)


def get_gene_records(args, gene_name):
    '''Retrieve the gene XML records for a given gene from the
    gene cache, or from online.
    '''
    if args.online:
        gene_ids = get_gene_ids(args, gene_name)
        for gene_id in gene_ids:
            gene_records_xml = fetch_records_from_ids(gene_ids)
            yield gene_id, gene_records_xml
    else:
        for gene_db_id, gene_file_path in \
            lookup_gene_cache_iter(args, gene_name):
            try:
                with open(gene_file_path) as gene_xml_file:
                    gene_xml_contents = gene_xml_file.read()
            except EnvironmentError as exception:
                exit("{}: failed to open gene XML file: {}".
                    format(gene_file_path, exception))
            yield gene_db_id, gene_xml_contents


def get_gene_ids(args, gene_name):
    '''Retrieve the NCBI gene database IDs for a gene name'''
    organism = args.organism
    search_term = '{}[organism] AND {}[sym]'.format(organism, gene_name)
    request = Entrez.esearch(db='gene', term=search_term, retmax=1000)
    try:
        result = Entrez.read(request)['IdList']
    except Exception as exception:
        logging.warn("Cannot get gene ids for gene {}: {}".
            format(gene_name, exception))
        return []
    else:
        return result

def fetch_records_from_ids(ids):
    '''Retrieve the NCBI gene database records for a list of database IDs'''
    database = 'gene'
    id_string = ','.join(list(ids))
    post_request = Entrez.epost(db=database, id=id_string)
    post_result = Entrez.read(post_request)
    web_env = post_result['WebEnv']
    query_key = post_result['QueryKey']
    fetch_request = Entrez.efetch(db=database,
                                 webenv=web_env,
                                 query_key=query_key,
                                 retmode='xml')
    return fetch_request.read()

class Regex(object):
    '''Regular expression search term'''

    def __init__(self, term):
        self.term = term
        self.match = re.compile(term)

    def search(self, string):
        """find the pattern in the string and
        return the context around the match"""
        spans = []
        for match in self.match.finditer(string):
            spans.append(match.span())
        if len(spans) > 0:
            return TermMatch(string, spans)
        else:
            return None

    def __str__(self):
        return self.term


class TermMatch(object):
    '''A representation of a term match containing its
    context (the whole string where the match occurred, plus
    all the spans (positions) of the matches'''
    def __init__(self, context, spans):
        self.context = context # string that was searched in
        self.spans = spans


def parse_search_term(string):
    '''Parse a string into a search term'''
    search_term = string.strip()
    return Regex(search_term)
