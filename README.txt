--------------------------------------------------------------------------------
Annokey - Annotation of gene lists with keyword hits from the NCBI gene
database and PubMed abstracts 
--------------------------------------------------------------------------------

Version: 1.0.0

Authors: Daniel Park (1), Sori Kang (3), Tú Nguyen-Dumont (1), Bernard J Pope (2,3)

         (1) Genetic Epidemiology Laboratory, Department of Pathology,
             The University of Melbourne.
         (2) Victorian Life Sciences Computation Initiative (VLSCI).
         (3) Department of Computing and Information Systems,
             The University of Melbourne.
         

Web page:       http://bjpop.github.io/annokey/

Repository:     https://github.com/bjpop/annokey


License: BSD

Requirements: Python 2.7, biopython 1.62, lxml 3.1, html 1.16

--------------------------------------------------------------------------------
General description
--------------------------------------------------------------------------------

Annokey is a command line tool for annotating gene lists with the results of a
keyword (or key phrase) search of the NCBI Gene database and linked PubMed
article abstracts. Its purpose is to help users prioritise genes by relevance
to a domain of interest, such as "breast cancer" or "DNA repair" etc. The user
steers the search by specifying a ranked list of keywords and terms that are
likely to be highly correlated with their domain of interest.

--------------------------------------------------------------------------------
Command line usage:
--------------------------------------------------------------------------------

usage: annokey [-h] [--version] [--online] [--cachesnapshot FILE]
               [--organism ORGANISM] [--email EMAIL_ADDRESS] [--genecache DIR]
               [--pubmedcache DIR] [--terms FILE] [--genes FILE]
               [--log FILENAME] [--delimiter {comma,tab}] [--report FILENAME]

Search NCBI for genes of interest, based on concept-keyword search.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --online              Search gene information from online (NCBI).
  --cachesnapshot FILE  Populate the gene cache from downloaded XML snapshot
                        of NCBI gene database
  --organism ORGANISM   Name of the organism to search
  --email EMAIL_ADDRESS
                        Your email address. This is required by NCBI for
                        online queries. You do not need to supply an email
                        address for cached queries.
  --genecache DIR       Save a cache of the downloaded results from NCBI gene
                        into this directory
  --pubmedcache DIR     Save a cache of the downloaded results from NCBI
                        pubmed into this directory
  --terms FILE          The tab separated file containing the search terms to
                        be searched.
  --genes FILE          The tab separated file containing the gene information
                        including name of the gene, one gene name per line.
  --log FILENAME        log progress in FILENAME, defaults to annokey_log.txt
  --delimiter {comma,tab}
                        Delimiter for gene file.
  --report FILENAME     Save a detailed search report as HTML page, defaults
                        to annokey_report.html
