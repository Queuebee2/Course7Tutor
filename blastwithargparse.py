import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description=f'This is the helpsection of {__file__} ', formatter_class=RawTextHelpFormatter)

parser.add_argument("-v", "--verbose",
                    help="verbosity level: turn on prints\n1 = a bit \n2 = a bit more",
                    type=int)

args = parser.parse_args()
if args.verbose:
    print("verbose is on")
    
