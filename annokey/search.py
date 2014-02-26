import sys
import csv
from Bio import Entrez
import logging
from process_xml import GeneParser 
from genecache import (lookup_gene_cache_iter,
    make_gene_cache_dirname, save_gene_cache)
from report import (DEFAULT_REPORT_FILE, init_report_page, report_hits, write_report)
from search_term import (parse_search_term)
from name import program_name


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
