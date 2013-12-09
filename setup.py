#!/usr/bin/env python

from distutils.core import setup

setup(
    name='annokey',
    version='0.1.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['annokey'],
    scripts=['annokey/genetoxml.py'],
    entry_points={
        'console_scripts': ['annokey = annokey.annokey:main']
    },
    url='https://github.com/bjpop/annokey',
    license='LICENSE.txt',
    description=(
        'Annokey: Annotation of gene lists with keyword hits '
        'from the NCBI gene database and PubMed abstracts.'),
    long_description=(
        'Annokey is a tool for annotating gene lists with the '
        'results of a keyword search of the NCBI gene database '
        'and linked Pubmed article information. Its purpose is '
        'to help users prioritise genetic mutations ordered by '
        'relevance to a domain of interest, such as "breast cancer"'
        'or "DNA repair pathways" etcetera. The user steers the '
        'search by specifying a ranked list of keywords and phrases '
        'that are likely to be highly correlated with their domain '
        'of interest.'),
    install_requires=[
        "biopython >= 1.62",
        "lxml >= 3.1"
    ],
)
