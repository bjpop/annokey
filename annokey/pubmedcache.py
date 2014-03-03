'''
Support for cacheing pubmed article information in the
filesystem.
'''

import os
from .hash import stable_string_hash


def make_pubmed_cache_dirname(cachedir, pubmed_id):
    '''Return the path of the directory for containing
    a particular pubmed entry. We hash the ID of the
    pubmed entry, and take modulus with 256 to spread
    the entries over 256 directories.
    '''
    hash_dir = stable_string_hash(pubmed_id) % 256
    return os.path.join(cachedir, str(hash_dir))


def save_pubmed_cache(cachedir, pubmed_id, pubmed_xml):
    '''Save an XML file for a pubmed entry in the appropriate
    place in the cache.
    '''
    pubmed_cache_dir = make_pubmed_cache_dirname(cachedir, pubmed_id)
    if not os.path.exists(pubmed_cache_dir):
        os.makedirs(pubmed_cache_dir)
    pubmed_cache_filename = os.path.join(pubmed_cache_dir, pubmed_id)
    with open(pubmed_cache_filename, 'w') as cache_file:
        cache_file.write(pubmed_xml)
