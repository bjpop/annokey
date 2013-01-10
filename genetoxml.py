#!/bin/env python

'''NCBI Gene Downlading Tool

This program downloads gene data from the ftp server and \
converting data to xml file by using linux.gene2xml tool.
This program checks the last modified time of the databse in the server, \
when it tries to download the gene database from the server. \
If the same version of database is in the local directory alreay, \
this program does not download the database from server. \
If you want to download the database anyway, please use the --force option. \
If you already have the database file that you need to look up and \
want to convert the file, please use the --convertonly option.

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


class Namespace:
    pass


def print_verbose(level, message):
    '''Print a message if the level is less than the verbose level.'''
    if level <= verboseLevel:
        print message
    if verboseLevel >= 3 and message.startswith('  [Error Log]'):
        raise


def flush_verbose(level, message):
    '''Flush a message if the level is less than the verbose level.'''
    if level <= verboseLevel:
        print message,
        sys.stdout.flush()


def convert_genetoxml(program, infile, outfile):
    '''Convert infile to xml file named outfile using linux.gene2xml.

    The path of the linux.gene2xml should be passed as an argument. \
    The program should be ether the full path of linux.gene2xml or \
    linux.gene2xml if path of linux.gene2xml is in PATH.
    '''
    # If infile is not ags format,
    # linux.gene2xml would complain about it.
    extension = infile.rfind('.ags')
    if extension == -1:
        print_verbose(0, '[Error] '
                         'The file to be converted is not ags file (%s).' %
                         infile)
        return

    class GeneToXmlError(Exception):
        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

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
            print_verbose(0, '[Error] '
                             'Could not execute the program linux.gene2xml, '
                             'perhaps it is not installed in your PATH? '
                             'Please add the program path to PATH or '
                             'use --gene2xmlpath option.')
        else:
            print_verbose(0, 'Undefined error occurs while '
                             'converting data to xml.')
            print_verbose(0, '  [Error Log] ' + str(e.errno))
    else:
        print_verbose(2, '  >> Converting %s to %s is successful.'
                         % (infile, outfile))

    return outfile


class Ftp(object):
    '''Ftp class for representation of ftplib.FTP

    This class supports connect(), disconnect() for FTP site, and \
    download(), get_last_modified_time(), and get_size()  for a file \
    on the FTP site.
    '''

    def __init__(self, ftpSite):
        self.site = ftpSite
        self.ftp = None

    def connect(self):
        # Connect FTP site and login.
        if self.ftp is not None:
            print_verbose(0, '[Error] ftp connection is already estabished.')
            return

        try:
            print_verbose(1, 'Starts connecting to server (%s).' % self.site)
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
            print_verbose(0, 'Undefined error occurs while '
                             'connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        return self.ftp

    def get_last_modified_time(self, filepath):
        # Get the last modified time of the file on the server.
        # The server gives the result code and the value.
        # e.g) '213 20120101051112'
        if self.ftp is None:
            print_verbose(0, '[Error] ftp servier is not connected.')
            return

        class MDTMError(Exception):
            def __init__(self, value):
                self.value = value

            def __str__(self):
                return repr(self.value)

        mdServer = None
        try:
            print_verbose(1, 'Starts getting the last modified time of data.')
            modifiedTime = self.ftp.sendcmd('MDTM ' + filepath)
            mdTimeResult = modifiedTime.split()
            if mdTimeResult[0] == '213':  # 213 means successful
                mdServer = int(mdTimeResult[1])
            else:
                raise MDTMError('Fail to get the last modified time of data '
                                'from server. [%s]' % modifiedTime)
            print_verbose(2, '  >> Getting the last modified time of data '
                             'is successful.')
            print_verbose(2, '  >> The last modified time of the data '
                             'at the server is %s' % str(mdServer))

        except (ftplib.all_errors, MDTMError) as e:
            print_verbose(0, 'Error occurs while connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        except KeyboardInterrupt:
            print_verbose(0, 'User interrupt. Finished program.')

        except Exception as e:
            print_verbose(0, 'Undefined error occurs while '
                             'connecting the ftp server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        return mdServer

    def get_size(self, filepath):
        '''Get the size of the file on the server.'''
        if self.ftp is None:
            print_verbose(0, '[Error] ftp servier is not connected.')
            return

        dataSize = None
        try:
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
            print_verbose(0, 'Undefined error occurs while '
                             'getting file size from server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        return dataSize

    def download(self, filepath, outfilename):
        '''Download the file from the server and save it named outfilename'''
        if self.ftp is None:
            print_verbose(0, '[Error] ftp servier is not connected.')
            return

        storefile = outfilename
        print_verbose(1, 'Starts downoloading data from server.')

        try:
            # Get data size for progressbar
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
            print_verbose(0, '\nError occurs while'
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
            print_verbose(0, '\nUndefined error occurs '
                             'downloading data from server.')
            print_verbose(0, '  [Error Log] ' + str(e))

        return storefile

    def disconnect(self):
        '''Quit the connection.'''
        if self.ftp is not None:
            try:
                self.ftp.quit()
            # ftplib exception is always occured when
            # try to quit while downloading a file.
            # In this program, it always happens when
            # KeyboardInterrupt is occured while downloading.
            # As this exception is not the exceptions that should be handled,
            # ignore the exceptions.
            except ftplib.all_errors as e:
                pass
                #print_verbose(3, 'Error occurs while '
                #                 'quiting ftp connection.')
                #print_verbose(3, '  [Error Log] ' + str(e))
            except KeyboardInterrupt:
                pass
                #print_verbose(3, 'User interrupt. Finished program.')
            except Exception as e:
                pass
                #print_verbose(3, 'Undefined error occurs while '
                #                 'quitting ftp connection.')
                #print_verbose(3, '  [Error Log] ' + str(e))
            finally:
                self.ftp = None


def parse_args():
    parser = ArgumentParser(
        description=('This program downloads gene data from '
                     'the ftp server and converting data to '
                     'xml file by using linux.gene2xml tool. '
                     'This program checks the last modified time '
                     'of the databse in the server, '
                     'when it tries to download the gene database '
                     'from the server. '
                     'If the same version of database is in '
                     'the local directory alreay, '
                     'this program does not download the database '
                     'from server. '
                     'If you want to download the database anyway, '
                     'please use the --force option. '
                     'If you already have the database file that '
                     'you need to look up and '
                     'want to convert the file, '
                     'please use the --convertonly option.'))

    parser.add_argument('--verbose', metavar='VERBOSE', type=int, default=2,
                        help='Prints progress messages. '
                             'The users have to pass the level of messages. '
                             '0: No progress messages. '
                             '1: Prints what stage the program is on. '
                             '2: Prints detail information on progress. '
                             '3: Prints callstack when error occurs.')

    parser.add_argument('--ftpsite', metavar='FTPSITE', type=str,
                        help=('The ftp site to be connected for '
                              'downloading gene database. '
                              'Default is NCBI server '
                              '(ftp.ncbi.nlm.nih.gov).'))

    parser.add_argument('--filepath', metavar='GENEFILEPATH', type=str,
                        help=('The full file path in the ftp site that '
                              'would be downloaded. '
                              'The path should start from the root of '
                              'the ftp site. '
                              'For example, '
                              'if you want to download gene database of '
                              'Archaea from NCBI FTP site, you should pass '
                              '/gene/DATA/ASN_BINARY/Archaea_Bacteria/'
                              'Archaea.ags.gz as an arugment value. '
                              'Default is Homo spiens '
                              '(/gene/DATA/ASN_BINARY/Mammalia/'
                              'Homo_sapiens.ags.gz)'))

    parser.add_argument('--force', action='store_true',
                        help='Downloads the latest version of data from '
                             'the server without regarding '
                             'the existing file on directory. '
                             'If there is the same version in '
                             'local directory, the file would be overwritten.')

    parser.add_argument('--downloadonly', action='store_true',
                        help='Downloads the latest version of data from '
                             'the server but not converts the file to xml.')

    parser.add_argument('--convertonly', metavar='GENEFILE', type=str,
                        help=('Converts the input file to xml.'))

    parser.add_argument('--gene2xmlpath', metavar='PROGRAMPATH', type=str,
                        help=('The path where the gene2xml is in. '
                              'If the path for linux.gene2xml is not in '
                              'your PATH, you should pass this argument. '
                              'For example, if linux.gene2xml is in '
                              '/usr/local/gene2xml/1.3/bin and not in PATH '
                              ', then you should pass either '
                              '/usr/local/gene2xml/1.3/bin/linux.gene2xml '
                              'or /usr/local/gene2xml/1.3/bin.'))

    return parser.parse_args()


def main():

    def get_program_path(path):
        '''Determine the linux.gene2xml path by looking up the input path.'''
        gene2xml = None
        if path is not None:
            if path.endswith('linux.gene2xml'):
                gene2xml = path
            else:
                if path.endswith('/'):
                    gene2xml = '%slinux.gene2xml' % path
                else:
                    gene2xml = '%s/linux.gene2xml' % path
            if not os.path.exists(gene2xml):
                return None
        else:
            gene2xml = 'linux.gene2xml'
        return gene2xml

    # Parse the optional arguments.
    args = parse_args()
    # Set the verboseLevel for printing messages according to its level.
    global verboseLevel
    verboseLevel = args.verbose
    # The options downloadonly and convertonly cannot coexist together.
    if args.downloadonly and args.convertonly:
        print_verbose(0, '[Error]'
                         'Either --downloadonly or --convertonly '
                         'can be used at once')
        return
    # If the gene2xmlpath is given
    # but linux.gene2xml is not in the given path, finish the program.
    gene2xml = get_program_path(args.gene2xmlpath)
    if gene2xml is None:
        print_verbose(0, '[Error] The program path %s does not exist.' %
                         args.gene2xmlpath)
        return
    # Determine the ftp server.
    if args.ftpsite is None:
        ftpServer = ncbiServer
    else:
        ftpServer = args.ftpsite
    # Determine data to be downloaded.
    if args.filepath is None:
        dataFilepath = homoSapiens
    else:
        dataFilepath = args.filepath

    doConvert = True
    isUptodateXml = False
    uptodateAgsfile = None
    uptodateXmlfile = None
    # If convertonly is on, not downloading file.
    # If convertonly is off, start connecting to ftp server for downloading.
    if not args.convertonly:
        # Connect to ftp server and login.
        ftp = Ftp(ftpServer)
        if ftp.connect() is None:
            return

        # Check the last modifited time of the file.
        isUptodateAgs = True
        mdServer = ftp.get_last_modified_time(dataFilepath)
        if mdServer is None:
            ftp.disconnect()
            return

        # Extract file name and extension from filepath.
        # As downloaded file has the name in
        # {gene database name}_yyyymmddhhmmss.{extension} format,
        # extracting file name and extension is necessary.
        filename = dataFilepath[dataFilepath.rfind('/')+1:]
        extPosition = filename.find('.')
        filePrefix = filename[:extPosition]
        fileExt = filename[extPosition:]
        uptodateAgsfile = '%s_%s%s' % (filePrefix, mdServer, fileExt)
        uptodateXmlfile = '%s_%s%s' % (filePrefix, mdServer, '.xml')
        # If option force is not on,
        # check whether the local directory has the same version or not.
        if not args.force:
            # Get modified time of the files that were retrieved.
            mdFiles = []
            files = glob.glob('./%s_[0-9]*%s' % (filePrefix, fileExt))
            for file in files:
                m = re.match('\.\/%s_([0-9]{14})%s' %
                             (filePrefix, fileExt), file)
                mdFiles.append(int(m.group(1)))

            # Compare the modified time.
            # If there is no files starting with the prefix of database or
            # the file version is not the latest version,
            # the local directory is regarded as
            # not having the latest version.
            mdFiles.sort()
            if not mdFiles or mdFiles[-1] < mdServer:
                isUptodateAgs = False

        # Check whether the up-to-date version of xml is
        # in the local directory
        isUptodateXml = os.path.exists(uptodateXmlfile)

        # If option force is on or
        # if the local directory have not the latest version of the databse,
        # download the database from the server.
        # Otherwise, check whether corresponding xml file exists or not.
        # If the xml file is in the directory,
        # disconnect ftp connection and finish the program.
        # If the xml file is not in the directory, convert the arg file.
        if not isUptodateAgs or args.force:
            uptodateAgsfile = ftp.download(dataFilepath, uptodateAgsfile)
            if uptodateAgsfile is None:
                ftp.disconnect()
                return
        else:
            print_verbose(0, 'You already have the latest version of '
                             'the file (%s).' % uptodateAgsfile)

        ftp.disconnect()

    # Determine converting operation.
    if args.convertonly:
        # Check the corresponding xml file of the requested ags file.
        uptodateAgsfile = args.convertonly
        uptodateXmlfile = '%s%s' % \
                          (uptodateAgsfile[:uptodateAgsfile.rfind('.ags')],
                           '.xml')
        if os.path.exists(uptodateXmlfile):
            doConvert = False
            print_verbose(0, 'You already have the converted file (%s).' %
                             uptodateXmlfile)
        else:
            doConvert = True
    else:
        # Check the corresponding xml file of the downloaded ags file.
        if args.downloadonly:
            doConvert = False
        else:
            if isUptodateXml:
                doConvert = False
                print_verbose(0, 'You already have the latest version '
                                 'of the xml file (%s).' % uptodateXmlfile)
            else:
                doConvert = True

    # Convert the file to xml using linux.gene2xml.
    # If downloadonly option is give, converting would not be executed.
    if not args.downloadonly and doConvert:
        uptodateXmlfile = convert_genetoxml(gene2xml,
                                            uptodateAgsfile,
                                            uptodateXmlfile)
        if uptodateXmlfile is None:
            return

    # Done. Finish the program.
    print_verbose(1, 'Done.')


if __name__ == '__main__':
    main()
