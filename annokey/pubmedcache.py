import os
from hash import stable_string_hash
from lxml import etree
from StringIO import StringIO

def make_pubmed_cache_dirname(cachedir, pubmed_id):
    hash_dir = stable_string_hash(pubmed_id) % 256
    return os.path.join(cachedir, str(hash_dir))

def save_pubmed_cache(cachedir, pubmed_id, pubmed_xml):
    pubmed_cache_dir = make_pubmed_cache_dirname(cachedir, pubmed_id)
    if not os.path.exists(pubmed_cache_dir):
        os.makedirs(pubmed_cache_dir)
    pubmed_cache_filename = os.path.join(pubmed_cache_dir, pubmed_id)
    with open(pubmed_cache_filename, 'w') as cache_file:
        cache_file.write(pubmed_xml)
