import os
import sys
from ruffus import *
import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as iotools

PARAMS = P.get_parameters([
    "%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "../pipeline.yml",
    "pipeline.yml"
    ])

# assembly {{{
@follows(
    mkdir("fastqs"),
    mkdir("contigs")
    )
#@transform(
#    "fastqs/*",
#    regex(".*"),
#    "contigs"
#    )
@split(
    "fastqs",
    "contigs/*.fa"
    )
def assembly(infile, outfile):
    ''' Contig assembly '''
    pass

# }}}

# fsm {{{
@merge(
    assembly,
    "kmers.gz"
    )
def fsm(infile, outfile):

    ''' Kmer mining/ counting with fsm-lite '''

    to_cluster = True

    statement = '''
    ls contigs | awk -F. '{print $1 "\t" $0}' > contigs_list.txt &&
    cd contigs &&
    fsm-lite -l ../contig_list.txt -s 6 -S 610 -v -t fsm_kmers | gzip -c > ../%(outfile)s
    '''

    P.run(statement)

# }}}

# prokka {{{
@follows(
    mkdir("prokka")
    )
@transform(
    assembly,
    regex("contigs/(.*)\.fa"),
    r"prokka/\1\.gff",
    r"\1"
    )
def prokka(infile, outfile, idd):
    
    ''' Annotate with prokka '''

    to_cluster = True

    statement = '''
    prokka --centre X --compliant %(infile)s --outdir prokka --force --prefix %(idd)s &&
    mv prokka/%(idd)s.gff prokka
    '''

    P.run(statement)

# }}}

# roary {{{
@follows(
    prokka,
    mkdir("roary")
    )
@merge(
    prokka,
    "roary/accessory_binary_genes.fa.newick"
    )
def roary(infile, outfile):

    ''' Make tree with Roary '''

    to_cluster = True

    statement = '''
    roary -f %(outfile)s -e -n -v -r %(infile)s/*.gff 
    '''

    P.run(statement)

    pass

# }}}

# distanceFromTree {{{
@transform(
    roary,
    regex(".*"),
    "distances.tsv"
    )
def distanceFromTree(infile, outfile):
    
    ''' Get distances from a phylogeny tree that has been midpoint rooted '''

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    statement = '''
    python %(PY_SRC_PATH)s/phylogeny_distance.py 
        --calc-C %(infile)s 
        > %(outfile)s
    '''
    
    P.run(statement)

# }}}

# splitPhenos {{{
@follows(
    mkdir("phenos")
    )
@split(
    "phenos.tsv",
    "phenos/*.tsv"
    )
def splitPhenos(infile, outfiles):
    ''' split a tsv file into multiple tsv files by column '''

    statement = '''
    cols=`awk -F"\\t" '{print NF; exit}' %(infile)s` &&
    for col in $(seq 2 $cols); do
        pheno=`awk -F"\\t" -v col=$col 'NR==1{print $col}' %(infile)s` &&
        awk -F"\\t" -v col=$col '{print $1"\\t"$col}' %(infile)s > phenos/${pheno}.tsv;
    done
    '''

    P.run(statement)

# }}}

# pyseer {{{
@follows(
    mkdir("pyseer")
    )
@transform(
    splitPhenos,
    regex("phenos/(.*)\.tsv"),
    add_inputs(distanceFromTree, fsm),
    r"pyseer/\1.assoc"
    )
def pyseer(infiles, outfile):

    pheno = infiles[0]
    distances = infiles[1]
    kmers = infiles[2]

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

    P.run(statement)

# }}}

# makeRefList {{{
@follows (
    mkdir("refs")
    )
@merge(
    [prokka, "refs/*"],
    regex("([^/]+)/(.*).gff"),
    "ref.txt",
    )
def makeRefList(infiles, outfile):

    ''' Make a list of references for kmer mapping '''

    statement = '''
    ls prokka | grep ".gff" | awk '{print $1"\t%(prokka_dir)/"$1"\tref" }' > %(outfile)s &&
    ls %(ref_dir)s  | awk '{print $1"\t%(ref_dir)/"$1"\tref" }' 
    '''

    P.run(statement)

# }}}

# mapKmers {{{
#@transform(
#    pyseer,
#    regex("pyseer.dir/(.*).assoc"),
#    r"map.dir/\1\.txt",
#    add_inputs = [makeRefList]
#    )
#def mapKmers(infiles, outfile):
#
#    infile = infiles[0][0]
#    ref_genomes = infiles[1][0]
#
#    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))
#
#    statement = '''
#    python %(PY_SRC_PATH)/summarise_annotations.py
#    '''
## }}}

# full {{{
@follows (
    pyseer
    )
def full():
    pass

# }}}

def main():
    P.main(sys.argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
