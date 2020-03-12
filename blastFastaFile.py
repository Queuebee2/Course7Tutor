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
# 


import argparse
from argparse import RawTextHelpFormatter

from Bio.Blast.NCBIWWW import qblast as BLAST
from Bio.Blast.NCBIXML import parse as parseXML
from Bio import Entrez


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
        print(f"trying to BLAST \n{kwargs['sequence'][:20]}...")
    
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
                if verbose:
                    print('****Alignment****')
                    print('sequence:', alignment.title)
                    print('length:', alignment.length)
                    print('score:', hsp.score)
                    print('gaps:', hsp.gaps)
                    print('e-value:', hsp.expect)
                    print(hsp.query[0:90] +'...')
                    print(hsp.match[0:90] +'...')
                    print(hsp.sbjct[0:90] +'...')
                    


# digest user input
args = parser.parse_args(["blastn","nt","ORFFOUND2.fa","milain.lambers@gmail.com","-o","blast_fromcommandline","-v","-e","1"])

# set email
Entrez.email = args.email_address

if args.verbose:
    print("verbose is on")
    print('these are the passed arguments:')
    print(args)


options = {'program': args.program,
           'database': args.database,
           'sequence': open(args.query).read()}

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
parseBLAST(xml_file_name=outputfile, verbose=args.verbose)

    


