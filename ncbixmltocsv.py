#!/bin/env python 
# this is a test comment.
# this is another test comment.

import lxml.etree as ET
import sys
import csv
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(
        description = 'Search NCBI for genes of interest, based on concept-keyword search')
    parser.add_argument(
        '--xml', metavar = 'INFILE', type=file,
        required=True, help ='name of the NCBI XML file')
    parser.add_argument(
        '--csv', metavar = 'OUTFILE', type=argparse.FileType('wb'),
        default=sys.stdout,
        help ='name of the output CSV file, defaults to stdout')
    return parser.parse_args()

def process_XML(inXML, writer):
    name = ''
    rifs = []
    commentary_type = None
    context = ET.iterparse(inXML)
    for action, elem in context:
        tag = elem.tag
        if tag == 'Gene-ref_locus':
            name = elem.text
        # Assumes commentary_type always appears before commentary_text
        elif tag == 'Gene-commentary_type':
            commentary_type = elem.get('value')
        elif tag == 'Gene-commentary_text' and commentary_type == 'generif':
            rifs.append('<' + elem.text + '>')
        elif tag == 'Entrezgene':
            tabbed_rifs = ('\t'.join(rifs))
            writer.writerow([name, tabbed_rifs])
            name = ''
            rifs = []
        elem.clear()

def main():
    args = parseArgs()
    writer = csv.writer(args.csv, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    process_XML(args.xml, writer)

if __name__ == '__main__':
    main()
