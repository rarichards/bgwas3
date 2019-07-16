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

PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))
R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

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
# getKmers {{{
@follows (
    getContigs
    )
@transform(
    input = getContigs,
    filter = regex("contigs"),
    output = r"kmers.gz"
    )
def getKmers(infile, outfile):
    ''' Kmer mining/ counting with fsm-lite '''

    statement = '''
    ls %(infile) | awk -F. '{print $1 "\t" $0}' > contig_list.txt &&
    fsm-lite
        -l contig_list.txt
        -s 6
        -S 610
        -v
        -t fsm_kmers
    | gzip 
        -c
        > %(outfile)
    '''

    P.run(statement)

# }}}
# getPhylogeny {{{
@follows(getContigs)
@transform(
    input = getContigs,
    filter = regex("contigs"),
    output = r"phylogeny.tree"
    )
def getPhylogeny(infile, outfile):
    
    ''' Get phylogeny '''

    pass

# }}}
# getDistances {{{
@follows(getPhylogeny)
@transform(
    input = getPhylogeny,
    filter = regex("phylogeny.tree"),
    output = r"distances.tsv"
    )
def getDistances(infile, outfile):
    
    '''Get distances from a phylogeny tree that has been midpoint rooted'''

    statement = '''
    python %(PY_SRC_PATH)s/phylogeny_distance.py 
        --calc-C %(infile)s 
        > %(outfile)s
    '''
    
    P.run(statement)

# }}}
# getPhenos {{{
@split(
    input = "pheno.tsv",
    output = r"phenos/.*\.tsv"
    )
def getPhenos(infile, outfile):

    statement = '''
    Rscript %(R_SRC_PATH)s/split_tsv.R %(infile)s phenos
    '''

    P.run(statement)

# }}}
# getAssoc {{{
@follows(
    getKmers,
    getDistances,
    getPhenos
    )
@transform(
    input = getPhenos,
    filter = regex("phenos/(.*)\.tsv"),
    output = r"\1\.assoc",
    add_inputs = [getDistances, getKmers]
    )
def get_assoc(infiles, outfile):

    phenos = infiles[0]
    distances = infiles[1][0]
    kmers = infiles[1][1]
    
    statement = '''
    pyseer
        --phenotypes=%(pheno)s
        --kmers=%(kmers)s
        --distances=%(distances)s
        --min-af=0.01
        --max-af=0.99
        --cpu=15
        --filter-pvalue=1E-8
        > %(outfile)s
    '''

# }}}
# getKmerAnnotation {{{
@follows(
    getAssoc,
    )
@transform(
    input = getAssoc,
    filter = regex("associations/(.*).assoc"),
    output = r"annotated/\1\.txt",
    add_inputs = ["ref_genomes"]
    )
def getKmerAnnotation(infiles, outfile):

    ref_genomes = infiles[1][0]

    statement = '''
    ls %(ref_genomes) | awk '{print $1 "\t" ref}' > ref_genome_list.txt &&
    python %(PY_SRC_PATH)/summarise_annotations.py
    '''
# }}}
# full {{{
@follows (
    getKmerAnnotation
    )
def full():
    pass

# }}}

def main():
    P.main(sys.argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
