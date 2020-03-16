# Date      6-3-2020
# Author    Milain Lambers
# GitHub    Queuebee2
############################################################
#
# inspiration https://github.com/Gurdhhu/bioinf_scripts
#
# Purpose
# commandline tool (in the making) to perform BLAST searches
# through the biopython connection to the NCBI blast api
#
# later on, some tricks are needed to save not as XML but
# tabular? or something? To make annotating stuff easier
# idk yet
#
# should probably take a look here
# https://github.com/peterjc/galaxy_blast/blob/master/tools/ncbi_blast_plus/blastxml_to_tabular.py
# 
#
# TODO
# - [ ] check validity of fastafile > throw error
# - [ ] check usability before storing (eval, type of prot)
# - [x] check existence in storage before blasting > else skip
# - [x] store more information ( which seq, what header.. )
#
# 
#
#
#
# 
#
# 
#
# 
#
#



import argparse
from argparse import RawTextHelpFormatter

from Bio.Blast.NCBIWWW import qblast as BLAST
from Bio.Blast.NCBIXML import parse as parseXML
from Bio import Entrez

DEFAULT_STORE_BLAST_CSV= "blast_results.csv"
DEFAULT_BLAST_OUTPUTNAME = "blast_output.xml"

parser = argparse.ArgumentParser(description=f"This is the helpsection of {__file__}",
                                 formatter_class=RawTextHelpFormatter)


# accepted arguments
parser.add_argument("-v", "--verbose",
                    help="""verbosity level: turn on prints""", action='store_true')

parser.add_argument('program',
                    help="""blast program to use.\n\ndefault = blastn\n\navailable options:\n\t- blastn\n\t- blastx""",
                    default='blastn',
                    choices=['blastn','blastx', 'blastp'])

parser.add_argument("database",
                    help="""database to query.\n\ndefault = nt\n\navailable options:\n\t- nt\n""",
                    default='nt',
                    choices=['nt'])

parser.add_argument("query",
                    help="""fastafile to query blast""")

parser.add_argument("-o", "-out", "--outputfile",
                    help="""the name of the output file""",
                    type=str)

parser.add_argument("-e", "--evalue",
                    type=int,
                    help="""the minimal e-value given as int""")

parser.add_argument("-ext", "--extention",
                    help="""disable use of .xml fileextention""", action='store_true')

parser.add_argument("email_address",
                    help="""your email address, this seems to be essential for the program to work (?!)""")

# parse added arguments





# helper functions

def doBlast(verbose=False, **kwargs):
    """ with a kwarg dict, blast a query sequence/sequences
        against the selected databse with a selected blast program
        and other passed parameters through the Bio.Blast.NCBIWWW.qblast
        method
    """

    if verbose:
        print('blast settings :')
        for k, v in kwargs.items():
            if k != 'sequence':
                print(f'{k} : {v}')
        print(30*"-")
        print(f"trying to BLAST \n{kwargs['sequence'][:25]}...")
    
    blast_result_handle = BLAST(**kwargs)
    
    if verbose:
        print(f"successfully did a BLAST with {kwargs['sequence'][:20]}")
        
    # execute blast return_handle
    return blast_result_handle
                                   

def saveXML_BLAST(blast_result_handle,filename=DEFAULT_BLAST_OUTPUTNAME, verbose=False):
    # open output file
    blast_file = open(filename, 'w')

    if verbose:
        print(f"saving results in {filename}")
    # write results ( xml-string )
    blast_file.write(blast_result_handle.read())
    blast_file.seek(0)
    blast_file.close()


def saveTABULAR_BLAST():
    """ placeholder  for
    implementing https://github.com/Gurdhhu/bioinf_scripts

    """

def parseBLAST(xml_file_name=DEFAULT_BLAST_OUTPUTNAME, verbose=False):
    """ THIS IS A GENERATOR!", yields blast results

    IMPORTANT TODO: set max amount of yields or based on evalue (or both )
    """
    # TODO patch verbose
    #
    # take a filename, asuming it is a textfile with xml formatted blast
    # results and output some standard information

    # NCBIXML.parse method can only be used on a file handle, strings dont work
    blast_file = open(xml_file_name, 'r')
    blast_records = parseXML(blast_file)

    
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                yield [ alignment.title,
                        alignment.length,
                        hsp.score,
                        hsp.expect]


def storeBLAST(blast_result, storagefile=DEFAULT_STORE_BLAST_CSV):
    with open(storagefile, 'a') as f:
        line = "\t".join([str(i) for i in blast_result])
        
        f. write(line+"\n")
        print(f' stored {line} ')
    
                    

def fastaReader(fastaFile):
    """Generator helper function to
    read and spit header-seq pairs one by one when iterated over"""

    print(f'opening {fastaFile} to read fastas..')
    with open(fastaFile, 'r') as f:

            header = f.readline()
            seq = ''

            # generator loop
            for line in f:
                    if line.startswith(">"):
                            yield header, seq
                            header = line
                            seq = ''
                    else:
                        seq += line

            # yield last pair
            yield header, seq



def exists(identifier, col=0, delim='\t', threshold=1, storagefile=DEFAULT_STORE_BLAST_CSV):
    """ identifier: an identifier for a sequence (e.g. a header)
        threshold : the amount of times the identifier can appear
        in the database ( a csv file for now )
        doesn't check whether lines are unique (as of now)

        the threshold was formerly based on the amount of blasts
        done per header, but this version is intended for
        one-blast-per sequence results
    """
    # this function can later be upgraded to query a database with
    # sql

    try:
        with open(storagefile, 'r') as savefile:
            count = 0
            for line in savefile:
                if line.split(delim)[col] == identifier:
                    count += 1
                if count >= threshold:
                    return True
                
    except FileNotFoundError:
        print("creating savefile")
        with open(storagefile, 'w') as f:
            f.write("header\ttitle\tlength\tscore\texpect\n")

        return exists(identifier)
    
    # if the item exists less than 'threshold' times, return False
    return False

    
# digest user input
# leave parser.parse_args() empty if ran from commandline to parse commandline args
# otherwise pass a list of args like this
# args = parser.parse_args(["blastn","nt","fastaFileToBlast","your@email.adress","-o","blastxmlFile","-v","-e","1"])
args = parser.parse_args()

# set email
Entrez.email = args.email_address

if args.verbose:
    print("verbose is on")
    print('these are the passed arguments:')
    print(args)


fastaSpitter = fastaReader(args.query)
# TODO build a test to test if the fasta is right before actually starting


# main BLAST loop, takes ages for big files.
for header, sequence in fastaSpitter:

    if exists(header.strip("\n")):
        # skip headers we've already blasted
        print(f'skipping {header[:35]}')
        continue
    
    options = {'program': args.program,
               'database': args.database,
               'sequence': header+sequence}

    if args.evalue:
        options['expect'] = args.evalue

    blast_result_handle = doBlast(verbose=args.verbose, **options)


    if args.outputfile:
        
        if not ".xml" in args.outputfile and args.extention:
            outputfile = args.outputfile+".xml"
        else:
            outputfile = args.outputfile

        print(f"saving BLSAT output as {outputfile}")
        
        
    else:
        print(f"no outputfile specified...! blast results will be saved in XML file: {DEFAULT_BLAST_OUTPUTNAME}")
        outputfile = DEFAULT_BLAST_OUTPUTNAME
        
        
    saveXML_BLAST(blast_result_handle, filename=outputfile, verbose=args.verbose)


    print("parsing BLAST for ya now")
    for blast_result in parseBLAST(xml_file_name=outputfile, verbose=args.verbose):
        storeBLAST([header.strip("\n")]  + blast_result)

    


