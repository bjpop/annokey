'''
A stable hash function for strings.
The built-in hash function in Python is not guaranteed to
produce the same results over different versions of Python.
We use the hash function to generate file paths, so we need
a version which is stable across versions of Python.
'''

import hashlib

def stable_string_hash(string):
    '''A hash function for strings based on MD5 hashing which should be stable
    across Python implementations, unlike the built-in hash function.
    It doesn't matter if this is slow because we don't call it often.
    '''
    return int(hashlib.md5(string).hexdigest(), 16)
