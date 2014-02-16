import hashlib
import os

def stable_string_hash(str):
    '''A hash function for strings based on MD5 hashing which should be stable
    across Python implementations, unlike the built-in hash function. It doesn't
    matter if this is slow because we don't call it often.'''
    return int(hashlib.md5(str).hexdigest(), 16)

def make_pubmed_cache_dirname(cachedir, pubmed_id):
    hash_dir = stable_string_hash(pubmed_id) % 256
    return os.path.join(cachedir, str(hash_dir))
