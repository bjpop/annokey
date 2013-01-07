#!/bin/env python


import ftplib
import glob
import re
import subprocess
import os
import sys

from argparse import ArgumentParser

'''
    This script file downloads Homo sapiens gene information from NCBI ftp server, and convert the retrieved binary file to xml file.
'''

# global variables
ftpServer = 'ftp.ncbi.nlm.nih.gov'
homo_sapiens = '/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz'


class MDTMError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class GeneToXmlError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Namespace: pass


def parse_args():
    parser = ArgumentParser( description=('Downloading gene data from ftp server and converting data to xml file.'))

    parser.add_argument(
                        '--verbose',
                        metavar='VERBOSE',
                        type=int,
                        default=60, 
                        help='Prints progress message. User have to pass the level of message. 0: No progress message. 30: Prints on which stage the program is on. 60: Prints detail information on progress. 100: Prints callstack when error occurs.')

    parser.add_argument(
                        '--force',
                        action='store_true',
                        help='Downloads the latest version of data from server without regarding the existing file on directory.')

    parser.add_argument(
                        '--species',
                        metavar='SPECIES',
                        type=str, 
                        help=('The file path of ftp server that need to be downloaded. Default is Homo spiens (/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz)'))

    parser.add_argument(
                        '--downloadonly',
                        action='store_true',
                        help='Downloads the latest version of data from server but not converts the file to xml.')

    parser.add_argument(
                        '--convertonly',
                        metavar='GENEFILE',
                        type=str,
                        help=('Converts the input file to xml.'))

    return parser.parse_args()


def main():

    def print_verbose(level, message):
        if level <= args.verbose:
            print message
        if args.verbose == 100 and message.startswith('  [Error Log]'):
            raise

    def flush_verbose(level, message):
        if level <= args.verbose:
            print message,
            sys.stdout.flush()

    # parse the optional arguments 
    args = parse_args()
   
    if args.downloadonly and args.convertonly:
        print_verbose(0, 'Either --downloadonly or --convertonly can be used at once')
        return

    # determine data to be download
    if args.species is None:
        dataFilePath = homo_sapiens
    else:
        dataFilePath = args.species
    
    dataFile = dataFilePath[dataFilePath.rfind('/')+1:]
    extPosition = dataFile.find('.')
    filePrefix = dataFile[:extPosition]
    fileExt = dataFile[extPosition:]

    if args.convertonly is not None:
        storeFile = args.convertonly
    else:
        ftp = None
        try:
            print_verbose(30, 'Starts conneting to server (%s).' % ftpServer)
            ftp = ftplib.FTP(ftpServer)
            ftp.login()
            print_verbose(60, '  >> Connecting to server is successful.')

            # get the last modified time of server data
            # the return value would be '213 20120101051112'
    
            print_verbose(30, 'Starts getting the last modified time of data.')
            modifiedTime = ftp.sendcmd('MDTM ' + dataFilePath) 
            mdTimeResult = modifiedTime.split()
            if mdTimeResult[0] == '213': # means successful
                mdServer = int(mdTimeResult[1])
            else: 
                raise MDTMError('Fail to get the last modified time of data'+ 
		    	    ' from server. [%s]' % modifiedTime)

            print_verbose(60, '  >> Getting the last modified time of data is successful.')
            print_verbose(60, '  >> The last modified time of the data at the server is %s' % str(mdServer))

        except (ftplib.all_errors, MDTMError) as e:
            if ftp is not None: ftp.quit()
            print_verbose(0, 'Error occurs while connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(e))
            return

        except KeyboardInterrupt:
            if ftp is not None: ftp.quit()
            print_verbose(0, 'User interrupt. Finished program.') 
            return
        except Exception as e:
            if ftp is not None: ftp.quit()
            print_verbose(0, 'Undefined error occurs while connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(e))
            return

        # get modified time of the files that were retrieved
 
        mdFiles = []
        files = glob.glob('./%s_[0-9]*%s' % (filePrefix, fileExt))
        for file in files:
            m = re.match('\.\/%s_([0-9]{14})%s' % (filePrefix, fileExt) , file)
            mdFiles.append(int(m.group(1)))


        # compare modified time
        # if there is no files starting with the prefix of database or
        # the file version is not the latest version, 
        # retrieve the latest version from the server.

        mdFiles.sort()

        if len(mdFiles) == 0 or mdFiles[-1] < mdServer or args.force:
            storeFile = './%s_%s%s' % (filePrefix, str(mdServer), fileExt)
            print_verbose(30, 'Starts downoloading data from server.')        
            # download data from the server, and save  
            try: 
                # get data size for progressbar
                ftp.sendcmd('TYPE i')
                dataSize = ftp.size(dataFilePath)
                dataSize = float(dataSize)
                print_verbose(60, '  >> The file size is %.2fM' % (dataSize/(1024**2)))
                ns = Namespace()
                ns.downloadedSize = 0
                with open(storeFile, 'wb') as f:
                    def callback(data):
                        f.write(data)
                        ns.downloadedSize += len(data) 
                        percent = ns.downloadedSize/dataSize*100
                        flush_verbose(60, '\b'*8 + '%6.2f%%' % percent)

                    flush_verbose(60, '  >> Downloading...               ') 
                    ftp.retrbinary('RETR %s' % dataFilePath, callback)
                

            except (IOError, ftplib.all_errors) as e:
                # print out error message and remove partial file
                if os.path.exists(storeFile): os.remove(storeFile)
                print_verbose(0, '\nError occurs while downloading data from server.') 
                print_verbose(0, '  [Error Log] ' + str(e))
                return
            except KeyboardInterrupt:
                print_verbose(0, '\nUser interrupt. Finished program.')
                if os.path.exists(storeFile): os.remove(storeFile)
                return
            except Exception as e:
                if os.path.exists(storeFile): os.remove(storeFile)
                print_verbose(0, '\nUndefined error occurs while downloading data from server.')
                print_verbose(0, '  [Error Log] ' + str(e))
                return
            finally: 
                ftp.quit()

            print_verbose(60, '  \n>> Downloading data is successful.')
            print_verbose(60, '  >> %s is stored.' % storeFile)
        
        else:
            ftp.quit()
            print_verbose(0, 'You already have the lastest version. %s' % str(mdFiles[-1]))
            print_verbose(30, 'Done.')
            return

    # convert file to xml
    # if --convertonly option is on, 
    # the output file name would be 'input file'.xml 
    if not args.downloadonly:

        if args.convertonly:
            fileName = storeFile[:storeFile.find('.ags.gz')]
            outFile = '%s.xgs' % fileName
        else: 
            outFile = './%s_%s.xgs' % (filePrefix, str(mdServer)) 
        print_verbose(30, 'Starts converting %s to %s.' % (storeFile, outFile))
      
        try:
            ret = subprocess.check_output(['linux.gene2xml', '-i', storeFile,
                     		           '-o', outFile, '-b', '-c'],
	 				 stderr=subprocess.STDOUT )
            if ret != '': # error occurs
                raise GeneToXmlError(ret)
            
        except (subprocess.CalledProcessError, GeneToXmlError) as e:
            # print out error message and remove partial files
            if os.path.exists(storeFile): os.remove(storeFile)
            if os.path.exists(outFile): os.remove(outFile)
            print_verbose(0, 'Error occurs while converting data.')
            print_verbose(0, '  [Error Log] ' + str(e)) 
       
        except KeyboardInterrupt:
            print_verbose(0, 'User interrupt. Finished program.')
            if os.path.exists(storeFile): os.remove(storeFile)
            if os.path.exists(outFile): os.remove(outFile)
        except Exception as e:
            print_verbose(0, '\nUndefined error occurs while converting data to xml.')
            print_verbose(0, '  [Error Log] ' + str(e))
            if os.path.exists(storeFile): os.remove(storeFile)
            if os.path.exists(outFile): os.remove(outFile)
        else:
            print_verbose(60, '  >> Converting %s to %s is successful.' % (storeFile, outFile)) 
    
      
    print_verbose(30, 'Done.')



if __name__ == '__main__':
    main()

