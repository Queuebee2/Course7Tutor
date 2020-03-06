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
from Bio import Entrez
parser = argparse.ArgumentParser(description=f"This is the helpsection of {__file__}",
                                 formatter_class=RawTextHelpFormatter)


# accepted arguments
parser.add_argument("-v", "--verbose",
                    help="""verbosity level: turn on prints""", action='store_true')

parser.add_argument('program',
                    help="""blast program to use.\n\ndefault = blastn\n\navailable options:\n\t- blastn\n\t- blastx""",
                    default='blastn',
                    choices=['blastn','blastx'])

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
                    help="""the minimal e-value """)

parser.add_argument("email",
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
        print(f"trying to BLAST \n{kwargs['sequence'][:20]}...")
    
    blast_result_handle = BLAST(**kwargs)
    
    if verbose:
        print(f"successfully did a BLAST with {kwargs['sequence'][:20]}")
        
    # execute blast return_handle
    return blast_result_handle
                                   

def saveXML_BLAST(blast_result_handle,filename="XML_BLAST_result", verbose=False):
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



# digest user input
args = parser.parse_args()

# set email
Entrez.email = args.email

if args.verbose:
    print("verbose is on")
    print('these are the passed arguments:')
    print(args)


options = {'program': args.program,
           'database': args.database,
           'sequence': open(args.query).read(),
           'matrix_name': "BLOSUM62", }

if args.evalue:
    options['expect'] = args.evalue

blast_result_handle = doBlast(verbose=args.verbose, **options)


if args.outputfile:
    saveXML_BLAST(blast_result_handle, filename="XML_BLAST_result")
    # open output file
    blast_file = open(filename, 'w')

    # write results ( xml-string )
    blast_file.write(blast_result_handle.read())
    blast_file.seek(0)
    blast_file.close()
else:
    print("no outputfile specified...! blast results lost in void!")



    


