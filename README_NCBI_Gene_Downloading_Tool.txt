NCBI Gene Downloading Tool

This program downloads gene data from the ftp server and 
converting data to xml file by using linux.gene2xml tool.
This program checks the last modified time of the database in the server,
when it tries to download the gene database from the server. 
If the same version of database is in the local directory already, 
this program does not download the database from server. 
If you want to download the database anyway, please use the --force option. 
If you already have the database file that you need to look up and 
want to convert the file, please use the --convertonly option.

The ftp site to be connected and the gene database to be downloaded are 
determined by the option arguments passed. The default values are
 - ftp site: NCBI ftp server (ftp.ncbi.nlm.nih.gov)
 - gene database: Homo sapiens 
                 (/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz)

Downloading gene file:
Before downloading data file from the server, it checks the last modified 
time of the database in the server. If the local directory (where the 
program is executed) has the same version, downloading does not executed. 
If the local directory has not the latest version, it downloads. However, 
by passing optional argument (--force), the database file would be downloaded
and converted without regarding the file in the local directory.

Converting gene file to xml:
This program uses linux.gene2xml tool as a converter. This program assume that
your PATH has the relevant path of linux.gene2xml. If your PATH does not have,
you should pass the option argument (--gene2xmlpath) to inform where 
linux.gene2xml is. For example, if linux.gene2xml is located at 
/usr/local/gene2xml/1.3/bin/linux.gene2xml and PATH does not have it, 
the argument should be 
--gene2xmlpath=/usr/local/gene2xml/1.3/bin/linux.gene2xml or 
--gene2xmlpath=/usr/local/gene2xml/1.3/bin


Decision Diagram 


                                   Start
                                     |
                Yes <---------- convert only?
                 |                   | No
                 |                   v 
                 |   Yes <---- force download?
                 |    |              | No
                 |    |              v
                 |    |  is the latest version in directory? -----> Yes
                 |    |              | No                          |
                 |    |              v                             |
                 |     ---->   Download database                   |
                 |                   |                             |
                 |                   | <---------------------------
                 |                   |              
                 |              download only?   ----------------> Yes
                 |                   | No                          |
                 |                   v                             |
                 -------->  is the xml in directory? ------> Yes   |
                                     | No                    |     |
                                     v                       |     |
                                Convert to xml               |     |
                                     |                       |     |
                                   Exit  <-------------------------     


Options
        [-h] [--verbose VERBOSE] [--ftpsite FTPSITE]
        [--filepath FILEPATH] [--force] [--downloadonly]
        [--convertonly GENEFILE] [--gene2xmlpath PROGRAMPATH]
	
        --verbose: Prints progress messages. The level of detail of the 
                   message could be adjustable by passing arguments. 
                   The levels are 0, 1, 2, and 3. The default level is 2.
                   Level 0: Prints out nothing except error message.
                   Level 1: Prints out the brief status of the progress.
                   Level 2: Prints out the detail of the progress.
                   Level 3: Same with Level 2 but displays callstack when 
                            error occurs.
         --------------------------------------------------------------
                error msg  | brief status | detail status | callstack
           0  |     o      |       x      |       x       |     x
           1  |     o      |       o      |       x       |     x
           2  |     o      |       o      |       o       |     x
           3  |     o      |       o      |       o       |     o
         --------------------------------------------------------------
         e.g) --verbose=2
         Starts connecting to server (ftp.ncbi.nlm.nih.gov).
           >> Connecting to server is successful.
         Starts getting the last modified time of data.
           >> Getting the last modified time of data is successful.
           >> The last modified time of the data at the server is 
              20130106185643
         Starts downloading data from server.
           >> The file size is 135M
           >> Downloading...         10.34% 
          
        --ftpsite: The ftp site to be connected for downloading gene database.
                   Default is NCBI server (ftp.ncbi.nlm.nih.gov).
         e.g) --ftpsite=ftp.ncbi.nlm.nih.gov

        --filepath: The full file path in the ftp site that would be 
                    downloaded. The path should start from the root of the 
                    ftp site. For example, if you want to download gene 
                    database of Archaea from NCBI FTP site, you should pass 
                    /gene/DATA/ASN_BINARY/Archaea_Bacteria/Archaea.ags.gz as 
                    an argument value. Default is Homo sapiens 
                    (/gene/DATA/ASN_BINARY/Mammalia/Homo_sapiens.ags.gz)
         e.g) --filepath=/gene/DATA/ASN_BINARY/Mammalia/Bos_taurus.ags.gz
       
        --force: Downloads database file without regarding the file in the 
                 local directory. If there is the same version in the local 
                 directory, the file would be overwritten. 
         e.g) --verbose=1 --force

        --downloadonly: Downloads the latest version of data from server but 
                        not converts the file to xml. 
         e.g) --verbose=60 --force --downloadonly
 
        --convertonly: Converts the input file to xml.
         e.g) --convertonly=Homo_sapiens_20130106185643.ags.gz

        --gene2xmlpath: The path where the gene2xml is in. If the path for 
                        linux.gene2xml is not in your PATH, you should pass 
                        this argument. For example, if linux.gene2xml is in 
                        /usr/local/gene2xml/1.3/bin and not in PATH, you
                        should pass either 
                        /usr/local/gene2xml/1.3/bin/linux.gene2xml or 
                        /usr/local/gene2xml/1.3/bin.
         e.g) --gene2xmlpath=/usr/local/gene2xml/1.3/bin


Output
	Two files should be stored when the program is finished successfully.
        If either --downloadonly or --convertonly is passed, only one of the
        two file would be stored.

        1. {gene database name}_{last modified time}.ags.gz 
        2. {gene database name}_{last modified time}.xml

        The format of last modified time is yyyymmddhhmmss. Note that the last
        modified time is the server’s time.

         e.g) For --ftpsite=ftp.ncbi.nlm.nih.gov 
                  --filepath=/gene/DATA/ASN_BINARY/Mammalia/Bos_taurus.ags.gz,
              you would get
              Homo_sapiens_20130106185643.ags.gz
              Homo_sapiens_20130106185643.xml


Performance 

        1. Time elapsed for downloading gene binary file: 660 secs 
           (varies depending on network/server status)
        2. Time elapsed for converting file: 120 secs

         • performed on merri
         • binary file size : 136M
         • xml file size: 8.3G


Remarks: Verified that the output xml file correctly works for ‘NCBI Gene Search Tool ’.


