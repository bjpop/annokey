#!/bin/env python

'''NCBI Gene Downloading Tool

This program downloads gene data from the ftp server and \
converts data to xml file by using linux.gene2xml tool.
This program checks the last modified time of the database in the server, \
when it tries to download the gene database from the server. \
If the same version of database is in the local directory already, \
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
NCBI_SERVER = 'ftp.ncbi.nlm.nih.gov'
HOMO_SAPIENS = '/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz'
VERBOSE_LEVEL = 2


class Namespace(object):
    '''Class introducing a new namespace'''
    pass


def print_verbose(level, message):
    '''Print a message if the level is less than the verbose level.'''
    if level <= VERBOSE_LEVEL:
        print message
    if VERBOSE_LEVEL >= 3 and message.startswith('  [Error Log]'):
        raise


def flush_verbose(level, message):
    '''Flush a message if the level is less than the verbose level.'''
    if level <= VERBOSE_LEVEL:
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
        '''Exception raised on conversion of gene entry to XML'''
        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    print_verbose(1, 'Started converting %s to %s.' % (infile, outfile))
    print_verbose(2, '  >> Converting ...')
    try:
        ret = subprocess.check_output([program, '-i', infile,
                                       '-o', outfile, '-b', '-c'],
                                      stderr=subprocess.STDOUT)
        if ret != '':  # error occurs
            raise GeneToXmlError(ret)

    except (subprocess.CalledProcessError, GeneToXmlError) as exception:
        if os.path.exists(outfile):
            os.remove(outfile)
        outfile = None
        print_verbose(0, 'Error occurred while converting data.')
        print_verbose(0, '  [Error Log] ' + str(exception))

    except KeyboardInterrupt:
        if os.path.exists(outfile):
            os.remove(outfile)
        outfile = None
        print_verbose(0, 'User interrupt. Finished program.')

    except Exception as exception:
        if os.path.exists(outfile):
            os.remove(outfile)
        outfile = None
        if exception.errno == 2:  # linux.gene2xml is not in PATH
            print_verbose(
                0,
                '[Error] '
                'Could not execute the program linux.gene2xml.'
                '\n        Perhaps it is not installed in your PATH? '
                '\n        Please add the program path to PATH or '
                'use --gene2xmlpath option.')
        else:
            print_verbose(0, 'Undefined error occurs while '
                             'converting data to xml.')
            print_verbose(0, '  [Error Log] ' + str(exception.errno))
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
        '''Connect FTP site and login.'''

        if self.ftp is not None:
            print_verbose(0, '[Error] ftp connection is already established.')
            return

        try:
            print_verbose(1, 'Starts connecting to server (%s).' % self.site)
            self.ftp = ftplib.FTP(self.site)
            self.ftp.login()
            print_verbose(2, '  >> Connecting to server is successful.')

        except ftplib.all_errors as exception:
            self.disconnect()
            print_verbose(0, 'Error occurs while connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(exception))

        except KeyboardInterrupt:
            self.disconnect()
            print_verbose(0, 'User interrupt. Finished program.')

        except Exception as exception:
            self.disconnect()
            print_verbose(0, 'Undefined error occurs while '
                             'connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(exception))

        return self.ftp

    def get_last_modified_time(self, filepath):
        '''Get the last modified time of the file on the server.
        The server gives the result code and the value.
        e.g) '213 20120101051112'
        '''
        if self.ftp is None:
            print_verbose(0, '[Error] ftp server is not connected.')
            return

        class MDTMError(Exception):
            '''Representation of modification time'''
            def __init__(self, value):
                self.value = value

            def __str__(self):
                return repr(self.value)

        md_server = None
        try:
            print_verbose(1, 'Starts getting the last modified time of data.')
            modified_time = self.ftp.sendcmd('MDTM ' + filepath)
            md_time_result = modified_time.split()
            if md_time_result[0] == '213':  # 213 means successful
                md_server = int(md_time_result[1])
            else:
                raise MDTMError('Fail to get the last modified time of data '
                                'from server. [%s]' % modified_time)
            print_verbose(2, '  >> Getting the last modified time of data '
                             'is successful.')
            print_verbose(2, '  >> The last modified time of the data '
                             'at the server is %s' % str(md_server))

        except (ftplib.all_errors, MDTMError) as exception:
            print_verbose(0, 'Error occurs while connecting ftp server.')
            print_verbose(0, '  [Error Log] ' + str(exception))

        except KeyboardInterrupt:
            print_verbose(0, 'User interrupt. Finished program.')

        except Exception as exception:
            print_verbose(0, 'Undefined error occurs while '
                             'connecting the ftp server.')
            print_verbose(0, '  [Error Log] ' + str(exception))

        return md_server

    def get_size(self, filepath):
        '''Get the size of the file on the server.'''
        if self.ftp is None:
            print_verbose(0, '[Error] ftp server is not connected.')
            return

        data_size = None
        try:
            self.ftp.sendcmd('TYPE i')
            data_size = self.ftp.size(filepath)
            data_size = float(data_size)
            print_verbose(2, '  >> The file size is %.2fM'
                             % (data_size/(1024**2)))
        except ftplib.all_errors as exception:
            print_verbose(0, 'Error occurs while getting file size.')
            print_verbose(0, '  [Error Log] ' + str(exception))

        except KeyboardInterrupt:
            print_verbose(0, 'User interrupt. Finished program.')

        except Exception as exception:
            print_verbose(0, 'Undefined error occurs while '
                             'getting file size from server.')
            print_verbose(0, '  [Error Log] ' + str(exception))

        return data_size

    def download(self, filepath, outfilename):
        '''Download the file from the server and save it named outfilename'''
        if self.ftp is None:
            print_verbose(0, '[Error] ftp servier is not connected.')
            return

        storefile = outfilename
        print_verbose(1, 'Starts downloading data from server.')

        try:
            # Get data size for progressbar
            self.ftp.sendcmd('TYPE i')
            data_size = self.ftp.size(filepath)
            data_size = float(data_size)
            print_verbose(2, '  >> The file size is %.2fM'
                             % (data_size/(1024**2)))
            namespace = Namespace()
            namespace.downloaded_size = 0
            with open(storefile, 'wb') as the_file:
                def callback(data):
                    '''write data and display progress'''
                    the_file.write(data)
                    namespace.downloaded_size += len(data)
                    percent = namespace.downloaded_size/data_size*100
                    flush_verbose(2, '\b'*8 + '%6.2f%%' % percent)

                flush_verbose(2, '  >> Downloading...               ')
                self.ftp.retrbinary('RETR %s' % filepath, callback)

            print_verbose(2, '\n  >> Downloading data is successful.')
            print_verbose(2, '  >> %s is stored.' % storefile)

        except (IOError, ftplib.all_errors) as exception:
            if os.path.exists(storefile):
                os.remove(storefile)
            storefile = None
            print_verbose(0, '\nError occurs while'
                             'downloading data from server.')
            print_verbose(0, '  [Error Log] ' + str(exception))

        except KeyboardInterrupt:
            if os.path.exists(storefile):
                os.remove(storefile)
            storefile = None
            print_verbose(0, '\nUser interrupt. Finished program.')

        except Exception as exception:
            if os.path.exists(storefile):
                os.remove(storefile)
            storefile = None
            print_verbose(0, '\nUndefined error occurs '
                             'downloading data from server.')
            print_verbose(0, '  [Error Log] ' + str(exception))

        return storefile

    def disconnect(self):
        '''Quit the connection.'''
        if self.ftp is not None:
            try:
                self.ftp.quit()
            # ftplib exception always occurs when
            # try to quit while downloading a file.
            # In this program, it always happens when
            # KeyboardInterrupt occurs while downloading.
            # As this exception is not the exceptions that should be handled,
            # ignore the exceptions.
            except ftplib.all_errors:
                pass
                #print_verbose(3, 'Error occurs while '
                #                 'quitting ftp connection.')
                #print_verbose(3, '  [Error Log] ' + str(e))
            except KeyboardInterrupt:
                pass
                #print_verbose(3, 'User interrupt. Finished program.')
            except Exception:
                pass
                #print_verbose(3, 'Undefined error occurs while '
                #                 'quitting ftp connection.')
                #print_verbose(3, '  [Error Log] ' + str(e))
            finally:
                self.ftp = None


def parse_args():
    '''Create a command line argument parser.'''

    parser = ArgumentParser(
        description=('This program downloads gene data from '
                     'the ftp server and converting data to '
                     'xml file by using linux.gene2xml tool. '
                     'This program checks the last modified time '
                     'of the database in the server, '
                     'when it tries to download the gene database '
                     'from the server. '
                     'If the same version of database is in '
                     'the local directory already, '
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
                              'Archaea.ags.gz as an argument value. '
                              'Default is Homo sapiens '
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
    '''Program entry point'''

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
    # Set the VERBOSE_LEVEL for printing messages according to its level.
    global VERBOSE_LEVEL
    VERBOSE_LEVEL = args.verbose
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
        print_verbose(0, '[Error] Could not find linux.gene2xml in %s' %
                         args.gene2xmlpath)
        return
    # Determine the ftp server.
    if args.ftpsite is None:
        ftp_server = NCBI_SERVER
    else:
        ftp_server = args.ftpsite
    # Determine data to be downloaded.
    if args.filepath is None:
        data_filepath = HOMO_SAPIENS
    else:
        data_filepath = args.filepath

    do_convert = True
    is_uptodate_xml = False
    uptodate_agsfile = None
    uptodate_xml_file = None
    # If convertonly is on, not downloading file.
    # If convertonly is off, start connecting to ftp server for downloading.
    if not args.convertonly:
        # Connect to ftp server and login.
        ftp = Ftp(ftp_server)
        if ftp.connect() is None:
            return

        # Check the last modified time of the file.
        is_uptodate_ags = True
        md_server = ftp.get_last_modified_time(data_filepath)
        if md_server is None:
            ftp.disconnect()
            return

        # Extract file name and extension from filepath.
        # As downloaded file has the name in
        # {gene database name}_yyyymmddhhmmss.{extension} format,
        # extracting file name and extension is necessary.
        filename = data_filepath[data_filepath.rfind('/')+1:]
        ext_position = filename.find('.')
        file_prefix = filename[:ext_position]
        file_ext = filename[ext_position:]
        uptodate_agsfile = '%s_%s%s' % (file_prefix, md_server, file_ext)
        uptodate_xml_file = '%s_%s%s' % (file_prefix, md_server, '.xml')
        # If option force is not on,
        # check whether the local directory has the same version or not.
        if not args.force:
            # Get modified time of the files that were retrieved.
            md_files = []
            files = glob.glob('./%s_[0-9]*%s' % (file_prefix, file_ext))
            for the_file in files:
                match = re.match(r'\.\/%s_([0-9]{14})%s' %
                             (file_prefix, file_ext), the_file)
                md_files.append(int(match.group(1)))

            # Compare the modified time.
            # If there is no files starting with the prefix of database or
            # the file version is not the latest version,
            # the local directory is regarded as
            # not having the latest version.
            md_files.sort()
            if not md_files or md_files[-1] < md_server:
                is_uptodate_ags = False

        # Check whether the up-to-date version of xml is
        # in the local directory
        is_uptodate_xml = os.path.exists(uptodate_xml_file)

        # If option force is on or
        # if the local directory have not the latest version of the database,
        # download the database from the server.
        # Otherwise, check whether corresponding xml file exists or not.
        # If the xml file is in the directory,
        # disconnect ftp connection and finish the program.
        # If the xml file is not in the directory, convert the arg file.
        if not is_uptodate_ags or args.force:
            uptodate_agsfile = ftp.download(data_filepath, uptodate_agsfile)
            if uptodate_agsfile is None:
                ftp.disconnect()
                return
        else:
            print_verbose(0, 'You already have the latest version of '
                             'the file (%s).' % uptodate_agsfile)

        ftp.disconnect()

    # Determine converting operation.
    if args.convertonly:
        # Check the corresponding xml file of the requested ags file.
        uptodate_agsfile = args.convertonly
        uptodate_xml_file = '%s%s' % \
                          (uptodate_agsfile[:uptodate_agsfile.rfind('.ags')],
                           '.xml')
        if os.path.exists(uptodate_xml_file):
            do_convert = False
            print_verbose(0, 'You already have the converted file (%s).' %
                             uptodate_xml_file)
        else:
            do_convert = True
    else:
        # Check the corresponding xml file of the downloaded ags file.
        if args.downloadonly:
            do_convert = False
        else:
            if is_uptodate_xml:
                do_convert = False
                print_verbose(0, 'You already have the latest version '
                                 'of the xml file (%s).' % uptodate_xml_file)
            else:
                do_convert = True

    # Convert the file to xml using linux.gene2xml.
    # If downloadonly option is give, converting would not be executed.
    if not args.downloadonly and do_convert:
        uptodate_xml_file = convert_genetoxml(gene2xml,
                                            uptodate_agsfile,
                                            uptodate_xml_file)
        if uptodate_xml_file is None:
            return

    # Done. Finish the program.
    print_verbose(1, 'Done.')


if __name__ == '__main__':
    main()
