import os
import sys
from ruffus import *
import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as iotools

P.get_parameters(
 ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
  "../pipeline.yml",
  "pipeline.yml"])

# getContigs {{{
@transform(
    input = "fastq",
    filter = regex("fastq"),
    output = r"contigs"
    )
def getContigs(infile, outfile):
    ''' Contig assembly '''
    pass

# }}}
# getContigsList {{{
@follows (
    getContigs
    )
@transform(
    input = getContigs,
    filter = regex("contigs"),
    output = r"list.txt"
    )
def getContigsList(infile, outfile):

    statement = '''
    ls %(infile) | awk -F. '{print $1 "\t" $0}' > %(outfile)
    '''

    P.run(statement)

# }}}
# full {{{
@follows(
    getContigsList
    )
def full():
    pass

# }}}

def main():
    P.main(sys.argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
