import os
import sys
from ruffus import *
import cgatcore.experiment as E
import cgatcore.iotools as iotools

def get_options(): # {{{
	import argparse

    parser = argparse.ArgumentParser(
			prog = "bgwas3" 
			description = "Description",
			)

    phenotypes = parser.add_argument_group('Test') # {{{

    phenotypes.add_argument(
			'--test',
            required=False,
            help='For testing (ignore)'
			)

	return parser.parse_args()

	# }}}

# }}}

def get_params(): # {{{
	from cgatcore import pipeline as P

	P.get_parameters(["pipeline.yml"])

	return P.PARAMS

# }}}

def main():
	options = get_options()
	params = get_params()

if __name__ == "__main__":
    main()
