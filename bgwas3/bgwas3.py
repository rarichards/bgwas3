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
    input = "fastq.dir",
    filter = regex("fastq.dir"),
    output = r"contigs.dir"
    )
def getContigs(infile, outfile):
    ''' Contig assembly '''
    pass

# }}}
# listContigs {{{
@follows (
    getContigs
    )
@transform(
    input = getContigs,
    filter = regex(".*"),
    output = r"contig_list.txt"
    )
def listContigs(infile, outfile):
    statement = '''
    ls %(infile)s | awk -F. '{print $1 "\t" $0}' > %(outfile)s
    '''
    P.run(statement)

# }}}
# getKmers {{{
@follows (
    getContigs
    )
@transform(
    input = getContigs,
    filter = regex(".*"),
    output = r"kmers.gz"
    )
def getKmers(infile, outfile):
    ''' Kmer mining/ counting with fsm-lite '''

    statement = '''
    ls %(infile)s | awk -F. '{print $1 "\t" $0}' > contig_list.txt &&
    cd %(infile)s &&
    ls &&
    fsm-lite -l ../contig_list.txt -s 6 -S 610 -v -t fsm_kmers | gzip -c > ../%(outfile)s
    '''

    P.run(statement)

# }}}
# prokka {{{
@follows(
    getContigs,
    mkdir("annotations.dir")
    )
@transform(
    "contigs.dir/*",
    regex("contigs.dir/(.*)\.fa"),
    r"annotations\.dir/\1/\1.gff",
    r"\1"
    )
def prokka(infile, outfile, idd):
    
    ''' Annotate with prokka '''

    to_cluster = True

    statement = '''
    prokka --centre X --compliant %(infile)s --outdir annotations.dir/%(idd)s --force
    '''

    P.run(statement)

# }}}
# roary {{{
@follows(
    prokka,
    )
@transform(
    input = "annotations.dir/*",
    filter = regex(r".*"),
    output = "roary.dir"
    )
def roary(infile, outfile):

    statement = '''
    roary -f %(outfile)s -e -n -v -r %(infiles)/*.gff 
    '''

    #P.run(statement)

    pass
# }}}

# getPhylogeny {{{
@follows(getContigs)
@transform(
    input = getContigs,
    filter = regex(".*"),
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

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    statement = '''
    python %(PY_SRC_PATH)s/phylogeny_distance.py 
        --calc-C %(infile)s 
        > %(outfile)s
    '''
    
    P.run(statement)

# }}}
# getPhenos {{{
@transform(
    input = "phenos.tsv",
    filter = regex("phenos.tsv"),
    output = r"phenos.dir"
    )
def getPhenos(infile, outfile):
    ''' split a tsv file into multiple tsv files by column '''

    statement = '''
    mkdir %(outfile)s &&
    cd %(outfile)s &&
    cols=`awk -F"\\t" '{print NF; exit}' ../%(infile)s` &&
    for col in $(seq 2 $cols); do
        pheno=`awk -F"\\t" -v col=$col 'NR==1{print $col}' ../%(infile)s` &&
        awk -F"\\t" -v col=$col '{print $1"\\t"$col}' ../%(infile)s > ${pheno}.tsv;
    done
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
def getAssoc(infiles, outfile):

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

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

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
