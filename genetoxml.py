#!/bin/env python

from ftplib import FTP

import glob
import re
import subprocess


'''
    This script file downloads Homo sapiens gene information from NCBI ftp server, and convert the retrieved binary file to xml file.
'''

homo_sapiens = '/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz'


def main():

    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()

    # get the last modified time of server data
    # the return value would be '213 20120101051112'
    
    modifiedTime = ftp.sendcmd('MDTM ' + homo_sapiens)
    mdTimeResult = modifiedTime.split()
    
    if mdTimeResult[0] == '213':  # means successful
        mdServer = int(mdTimeResult[1])
    else:  # retry??
        print 'Fail to get the last modified time.'
        ftp.quit()
        return

    # get modified time of the files that were retrieved
    
    mdFiles = []
    files = glob.glob('./Homo_sapiens_*.ags.gz')
    for file in files:
        m = re.match('\.\/Homo_sapiens_([0-9]*)\.ags\.gz', file)
        mdFiles.append(int(m.group(1)))

    # compare modified time
    # if there is no files starting with 'Homo_sapiens' or
    # the file version is not the latest version, 
    # retrieve the latest version from the server.
    
    mdFiles.sort()
    if len(mdFiles) == 0 or mdFiles[-1] < mdServer:
        storeFile = './Homo_sapiens_%s.ags.gz' % str(mdServer)
        f = open(storeFile, 'wb')

        def callback(data):
            f.write(data)
        
        ftp.retrbinary('RETR %s' % homo_sapiens, callback)
	f.close()
        
        # convert file to xml
        outFile = './Homo_sapiens_%s.xml' % str(mdServer)      
        subprocess.call(['linux.gene2xml', '-i', storeFile,
			 '-o', outFile, '-b', '-c'])
    else:
        print 'You have the lastest version. %s' % str(mdFiles[-1])


    ftp.quit()


if __name__ == '__main__':
    main()

