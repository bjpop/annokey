from lxml import etree
from StringIO import StringIO
import os
import hashlib
import logging
from process_xml import get_geneContent
from hash import stable_string_hash
from name import program_name
import sqlite3

def lookup_gene_cache_iter(args, gene_name):
    '''Gene records are stored in a file which is named using the
    NCBI gene databse ID.
    This generator searches for a gene name in the gene database and yields all
    the database entries that match that gene (there could be more than one).
    It returns the full path to the cached file plus the database entry ID.
    '''
    gene_cache_dir = make_gene_cache_dirname(args.genecache, args.organism, gene_name)
    try:
        dir_contents = os.listdir(gene_cache_dir)
    except OSError:
        logging.info("Could not find cache entry for {}".format(gene_name))
        return
    for geneDatabaseID in os.listdir(gene_cache_dir):
        filepath = os.path.join(gene_cache_dir, geneDatabaseID)
        yield geneDatabaseID, filepath 


def make_gene_cache_dirname(cachedir, organism, gene_name):
    # we normalise the gene name to upper case.
    gene_name_upper = gene_name.upper()
    organism_cache_dir = os.path.join(cachedir, organism)
    hash_dir = stable_string_hash(gene_name_upper) % 256
    return os.path.join(organism_cache_dir, str(hash_dir), gene_name_upper)

def save_gene_cache(args):

    #cachedir = args.genecache
    organism = args.organism

    conn = sqlite3.connect(args.genecache)
    cursor = conn.cursor()
    cursor.execute('''DROP TABLE IF EXISTS genes''')
    cursor.execute('''CREATE TABLE genes
                       (gene_id INTEGER PRIMARY KEY, name TEXT)''')

    try:
        with open(args.cachesnapshot) as xml_file:

            # read each "Entrezgene" record in the input XML and write it out to
            # a cache file. Sometimes one gene will have multiple entries. We store
            # each gene entry in a file based on its database ID.
            parser = etree.iterparse(xml_file, events=('end',), tag='Entrezgene')

            for event, elem in parser:
                gene_id, content = get_geneContent(elem)

                gene_name = ''
                gene_name_list = content['Gene Name']
                if len(gene_name_list) > 0:
                    gene_name = gene_name_list[0]
                if gene_name:
                    print('inserting {}'.format(gene_name))
                else:
                    print('########## NO NAME ##########')

                cursor.execute("INSERT INTO genes VALUES (?,?)",
                              (gene_id, gene_name))

                #gene_name_element = elem.find('.//Gene-ref_locus')

                #if gene_name_element is not None:
                #    gene_id, content = get_geneContent(elem)
                #    cursor.execute("INSERT INTO genes VALUES ({}, '{}')".format(gene_id, gene_name_element.text))
                #    print(gene_name_element.text)
                    #conn.commit()

                # find the official name of the gene
                #if gene_name_element is not None:
                #    gene_name = gene_name_element.text
                #    # find the database id of the gene
                #    gene_id = elem.find('.//Gene-track_geneid').text
                #    gene_cache_dir = make_gene_cache_dirname(cachedir, organism, gene_name)
                #    if not os.path.exists(gene_cache_dir):
                #        os.makedirs(gene_cache_dir)
                #    # Write a single 'Entrezgene' entry to a file in the cache using the
                #    # database ID for the file name
                #    gene_cache_filename = os.path.join(gene_cache_dir, gene_id)
                #    with open(gene_cache_filename, 'w') as cache_file:
                #        cache_file.write(etree.tostring(elem))

                # free up memory used by the XML iterative parser
                elem.clear()
                while elem.getprevious() is not None:
                   del elem.getparent()[0]

    except EnvironmentError as e:
        exit("{}: failed to open gene cache XML file: {}".format(program_name, e))
    finally:
        conn.commit()
        conn.close()
