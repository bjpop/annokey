#!/bin/env python

'''
    This script file downloads Homo sapiens gene information from NCBI ftp server, and convert the retrieved binary file to xml file.
'''

import os
import sys

import ftplib
import glob
import re
import subprocess

from argparse import ArgumentParser

# global variables
ncbiServer = 'ftp.ncbi.nlm.nih.gov'
homoSapiens = '/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz'
verboseLevel = 2


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


class Namespace:
    pass


def print_verbose(level, message):
    if level <= verboseLevel:
        print message
    if verboseLevel >= 3 and message.startswith('  [Error Log]'):
        raise


def flush_verbose(level, message):
    if level <= verboseLevel:
        print message,
        sys.stdout.flush()


def getProgramPath(arg):

    gene2xml = None
    if arg is not None:
        if arg.endswith('linux.gene2xml'):
            gene2xml = arg
        else:
            if arg.endswith('/'):
                gene2xml = '%slinux.gene2xml' % arg
            else:
                gene2xml = '%s/linux.gene2xml' % arg
        if not os.path.exists(gene2xml):
            return None
    else:
        gene2xml = 'linux.gene2xml'

    return gene2xml


# converts infile to xml file using linux.gene2xml.
# The path of the linux.gene2xml should be passed as an argument.
# It should be ether the full path or linux.gene2xml if PATH has its path.

def convertGenetoXml(program, infile):
    # the outfile has the same file name with the input file
    # but having xml extension.
    extension = infile.rfind('.ags')
    if extension == -1:  # file is not ags format
        print_verbose(0, '[Error]' +
                         'The file to be converted is not ags file (%s).'
                         % infile)
        return

    fileName = infile[:extension]
    outfile = '%s.xml' % fileName

    print_verbose(1, 'Starts converting %s to %s.' % (infile, outfile))
    try:
        ret = subprocess.check_output([program, '-i', infile,
                                       '-o', outfile, '-b', '-c'],
                                      stderr=subprocess.STDOUT)
        if ret != '':  # error occurs
            raise GeneToXmlError(ret)

    except (subprocess.CalledProcessError, GeneToXmlError) as e:
        if os.path.exists(outfile):
            os.remove(outfile)
        outfile = None
        print_verbose(0, 'Error occurs while converting data.')
        print_verbose(0, '  [Error Log] ' + str(e))

    except KeyboardInterrupt:
        if os.path.exists(outfile):
            os.remove(outfile)
        outfile = None
        print_verbose(0, 'User interrupt. Finished program.')

    except Exception as e:
        if os.path.exists(outfile):
            os.remove(outfile)
        outfile = None
        if e.errno == 2:  # linux.gene2xml is not in PATH
            print_verbose(0, '[Error]' +
                             'Could not execute the program linux.gene2xml,' +
                             'perhaps it is not installed in your PATH?' +
                             'Please add the program path to PATH or' +
                             'use --gene2xmlpath option.')
        else:
            print_verbose(0, 'Undefined error occurs' +
                             'while converting data to xml.')
            print_verbose(0, '  [Error Log] ' + str(e.errno))
    else:
        print_verbose(2, '  >> Converting %s to %s is successful.'
                         % (infile, outfile))
    finally:
        return outfile

# Ftp class supports
# connect(), disconnect(), download(), getLastModifiedTime(), and getSize().
# When operation finished, disconnect() should be called.


class Ftp(object):

    def __init__(self, ftpsite):
        self.site = ftpsite
        self.ftp = None

    def connect(self):
        if self.ftp is not None:
            print_verbose(0, '[Error] ftp connection is already estabished.')
            return

        try:
            print_verbose(1, 'Starts conneting to server (%s).' % self.site)
            self.ftp = ftplib.FTP(self.site)
            self.ftp.login()
            print_verbose(2, '  >> Connecting to server is successful.')

        except ftplib.all_errors as e:
            self.disconnect()
            print_verbose(0, 'Error occurs while connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        except KeyboardInterrupt:
            self.disconnect()
            print_verbose(0, 'User interrupt. Finished program.')

        except Exception as e:
            self.disconnect()
            print_verbose(0, 'Undefined error occurs' +
                             'while connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        finally:
            return self.ftp

    def getLastModifiedTime(self, filepath):
        # get the last modified time of server data
        # the return value would be '213 20120101051112'
        if self.ftp is None:
            print_verbose(0, '[Error] ftp servier is not connected.')
            return

        mdServer = None
        try:
            print_verbose(1, 'Starts getting the last modified time of data.')
            modifiedTime = self.ftp.sendcmd('MDTM ' + filepath)
            mdTimeResult = modifiedTime.split()
            if mdTimeResult[0] == '213':  # means successful
                mdServer = int(mdTimeResult[1])
            else:
                raise MDTMError('Fail to get the last modified time of data ' +
                                'from server. [%s]' % modifiedTime)
            print_verbose(2, '  >> Getting the last modified time of data ' +
                             'is successful.')
            print_verbose(2, '  >> The last modified time of the data ' +
                             'at the server is %s' % str(mdServer))

        except (ftplib.all_errors, MDTMError) as e:
            print_verbose(0, 'Error occurs while connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        except KeyboardInterrupt:
            print_verbose(0, 'User interrupt. Finished program.')

        except Exception as e:
            print_verbose(0, 'Undefined error occurs ' +
                             'while connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        finally:
            return mdServer

    def getSize(self, filepath):
        if self.ftp is None:
            print_verbose(0, '[Error] ftp servier is not connected.')
            return

        dataSize = None
        try:
            # get data size for progressbar
            self.ftp.sendcmd('TYPE i')
            dataSize = self.ftp.size(filepath)
            dataSize = float(dataSize)
            print_verbose(2, '  >> The file size is %.2fM'
                             % (dataSize/(1024**2)))
        except ftplib.all_errors as e:
            print_verbose(0, 'Error occurs while getting file size.')
            print_verbose(0, '  [Error Log] ' + str(e))

        except KeyboardInterrupt:
            print_verbose(0, 'User interrupt. Finished program.')

        except Exception as e:
            print_verbose(0, 'Undefined error occurs' +
                             'while downloading data from server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        finally:
            return dataSize

    def download(self, filepath, mdserver):
        # download file from the server and save
        if self.ftp is None:
            print_verbose(0, '[Error] ftp servier is not connected.')
            return

        # extract the file name and extension to make file name to be stored.
        filename = filepath[filepath.rfind('/')+1:]
        extPosition = filename.find('.')
        filePrefix = filename[:extPosition]
        fileExt = filename[extPosition:]
        storefile = './%s_%s%s' % (filePrefix, str(mdserver), fileExt)
        print_verbose(1, 'Starts downoloading data from server.')

        try:
            # get data size for progressbar
            self.ftp.sendcmd('TYPE i')
            dataSize = self.ftp.size(filepath)
            dataSize = float(dataSize)
            print_verbose(2, '  >> The file size is %.2fM'
                             % (dataSize/(1024**2)))
            ns = Namespace()
            ns.downloadedSize = 0
            with open(storefile, 'wb') as f:
                def callback(data):
                    f.write(data)
                    ns.downloadedSize += len(data)
                    percent = ns.downloadedSize/dataSize*100
                    flush_verbose(2, '\b'*8 + '%6.2f%%' % percent)

                flush_verbose(2, '  >> Downloading...               ')
                self.ftp.retrbinary('RETR %s' % filepath, callback)

            print_verbose(2, '  \n>> Downloading data is successful.')
            print_verbose(2, '  >> %s is stored.' % storefile)

        except (IOError, ftplib.all_errors) as e:
            if os.path.exists(storefile):
                os.remove(storefile)
            storefile = None
            print_verbose(0, '\nError occurs while' +
                             'downloading data from server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        except KeyboardInterrupt:
            if os.path.exists(storefile):
                os.remove(storefile)
            storefile = None
            print_verbose(0, '\nUser interrupt. Finished program.')

        except Exception as e:
            if os.path.exists(storefile):
                os.remove(storefile)
            storefile = None
            print_verbose(0, '\nUndefined error occurs' +
                             'while downloading data from server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        finally:
            return storefile

    def disconnect(self):
        if self.ftp is not None:
            self.ftp.quit()
            self.ftp = None


def parse_args():
    parser = ArgumentParser(description=('This program downloads gene data from the ftp server and converting data to xml file by using linux.gene2xml tool. This program checks the last modified time of the databse in the server, when it tries to download the gene database from the server. If the same version of database is in the local directory alreay, this program does not download the database from server. If you want to download the database anyway, please use the --force option. If you already have the database file that you need to look up and want to convert the file, please use the --convertonly option.'))

    parser.add_argument('--verbose', metavar='VERBOSE', type=int, default=2,
                        help='Prints progress messages. The users have to pass the level of messages. 0: No progress messages. 1: Prints what stage the program is on. 2: Prints detail information on progress. 3: Prints callstack when error occurs.')

    parser.add_argument('--ftpsite', metavar='FTPSITE', type=str,
                        help=('The ftp site to be connected for downloading gene database. Default is NCBI server (ftp.ncbi.nlm.nih.gov).'))

    parser.add_argument('--filepath', metavar='FILEPATH', type=str,
                        help=('The full file path in the ftp site that would be downloaded. The path should start from the root of the ftp site. For example, if you want to download gene database of Archaea from NCBI FTP site, you should pass /gene/DATA/ASN_BINARY/Archaea_Bacteria/Archaea.ags.gz as an arugment value. Default is Homo spiens (/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz)'))

    parser.add_argument('--force', action='store_true',
                        help='Downloads the latest version of data from server without regarding the existing file on directory. If there is the same version in local directory, the file would be overwritten.')

    parser.add_argument('--downloadonly', action='store_true',
                        help='Downloads the latest version of data from server but not converts the file to xml.')

    parser.add_argument('--convertonly', metavar='GENEFILE', type=str,
                        help=('Converts the input file to xml.'))

    parser.add_argument('--gene2xmlpath', metavar='PROGRAMPATH', type=str,
                        help=('The path where the gene2xml is in. If the path for linux.gene2xml is not in your PATH, you should pass this argument. For example, either /usr/local/gene2xml/1.3/bin/linux.gene2xml or /usr/local/gene2xml/1.3/bin would be good.'))

    return parser.parse_args()


def main():
    global verboseLevel

    # parse the optional arguments
    args = parse_args()
    verboseLevel = args.verbose

    if args.downloadonly and args.convertonly:
        print_verbose(0, '[Error]' +
                         'Either --downloadonly or --convertonly' +
                         'can be used at once')
        return
    # check the linux.gene2xml path
    gene2xml = getProgramPath(args.gene2xmlpath)
    if gene2xml is None:
        print_verbose(0, '[Error] The program path %s does not exist.' %
                         args.gene2xmlpath)
        return

    # determine the ftp server
    if args.ftpsite is None:
        ftpServer = ncbiServer
    else:
        ftpServer = args.ftpsite

    # determine data to be download
    if args.filepath is None:
        dataFilePath = homoSapiens
    else:
        dataFilePath = args.filepath

    storefile = None
    if args.convertonly is not None:
        storefile = args.convertonly
    else:
        # connect to ftp server and login
        ftp = Ftp(ftpServer)
        if ftp.connect() is None:
            return

        # check the last modifited time of the file,
        # if option --force is not give
        uptodatefile = True
        mdServer = ftp.getLastModifiedTime(dataFilePath)
        if mdServer is None:
            ftp.disconnect()
            return

        if not args.force:
            # check whether the local directory has the same version or not
            # get modified time of the files that were retrieved
            filename = dataFilePath[dataFilePath.rfind('/')+1:]
            extPosition = filename.find('.')
            filePrefix = filename[:extPosition]
            fileExt = filename[extPosition:]

            mdFiles = []
            files = glob.glob('./%s_[0-9]*%s' % (filePrefix, fileExt))
            for file in files:
                m = re.match('\.\/%s_([0-9]{14})%s' %
                             (filePrefix, fileExt), file)
                mdFiles.append(int(m.group(1)))

            # compare modified time
            # if there is no files starting with the prefix of database or
            # the file version is not the latest version,
            # retrieve the latest version from the server.
            mdFiles.sort()
            if len(mdFiles) == 0 or mdFiles[-1] < mdServer:
                uptodatefile = False

        # download the file
        if not uptodatefile or args.force:
            storefile = ftp.download(dataFilePath, mdServer)
            if storefile is None:
                ftp.disconnect()
                return
        else:
            print_verbose(0, 'You already have the latest version of' +
                             'the file (%s_%s%s).' %
                             (filePrefix, mdServer, fileExt))
            ftp.disconnect()
            return

        ftp.disconnect()
        # downloading is finished

    # convert the file to xml using linux.gene2xml
    # if --downloadonly option is give, converting would not be executed.
    if not args.downloadonly:
        outfile = convertGenetoXml(gene2xml, storefile)
        if outfile is None:
            return

    # done.
    # finish this program
    print_verbose(1, 'Done.')


if __name__ == '__main__':
    main()
