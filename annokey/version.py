'''
Export the version number of annokey.
The version number is obtained from the package.
'''

import pkg_resources  # part of setuptools

ANNOKEY_VERSION = pkg_resources.require("annokey")[0].version
