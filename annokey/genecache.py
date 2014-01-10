#!/bin/env python

from lxml import etree
from StringIO import StringIO
import os
import hashlib
import logging

def lookup_gene_cache_iter(args, gene_name):
    gene_cache_dir = make_gene_cache_dirname(args.genecache, args.organism, gene_name)
    try:
        dir_contents = os.listdir(gene_cache_dir)
    except OSError:
        logging.info("Could not find cache entry for {}".format(gene_name))
        return
    for filename in os.listdir(gene_cache_dir):
        file = '%s/%s' % (gene_cache_dir, filename)
        file = open(file)
        yield file
        file.close()

def stable_string_hash(str):
    '''A hash function for strings based on MD5 hashing which should be stable
    across Python implementations, unlike the built-in hash function. It doesn't
    matter if this is slow because we don't call it often.'''
    return int(hashlib.md5(str).hexdigest(), 16)

def make_gene_cache_dirname(cachedir, organism, gene_name):
    # we normalise the gene name to upper case.
    gene_name_upper = gene_name.upper()
    organism_cache_dir = os.path.join(cachedir, organism)
    hash_dir = stable_string_hash(gene_name_upper) % 256
    return os.path.join(organism_cache_dir, str(hash_dir), gene_name_upper)

def save_gene_cache(cachedir, organism, xml_file):

    # read each "Entrezgene" record in the input XML and write it out to
    # a cache file. Sometimes one gene will have multiple entries. We store
    # each gene entry in a file based on its database ID.
    #parser = etree.iterparse(StringIO(gene_records_xml), events=('end',), tag='Entrezgene')
    parser = etree.iterparse(xml_file, events=('end',), tag='Entrezgene')

    for event, elem in parser:
        # find the official name of the gene
        gene_name_element = elem.find('.//Gene-ref_locus')
        if gene_name_element is not None:
            gene_name = gene_name_element.text
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
        else:
            print("skipping item")
