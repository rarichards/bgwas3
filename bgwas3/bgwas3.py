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

# assembly {{{
@follows(
    mkdir("fastq.dir"),
    )
@transform(
    "fastq.dir",
    regex("fastq\.dir"),
    "contigs.dir"
    )
def assembly(infile, outfile):

    ''' Contig assembly '''

    pass
    #statement = '''
    #mkdir %(outfile)s
    #'''

    #P.run(statement)

# }}}
# fsm {{{
@follows (
    assembly
    )
@transform(
    "contigs.dir",
    regex("contigs\.dir"),
    "fsm.dir"
    )
def fsm(infile, outfile):

    ''' Kmer mining/ counting with fsm-lite '''

    to_cluster = True

    statement = '''
    mkdir fsm.dir &&
    ls %(infile)s | awk -F. '{print $1 "\t" $0}' > contig_list.txt &&
    cd %(infile)s &&
    fsm-lite -l ../contig_list.txt -s 6 -S 610 -v -t fsm_kmers | gzip -c > ../%(outfile)s/kmers.gz
    '''

    P.run(statement)

# }}}
# prokka {{{
@follows(
    assembly,
    mkdir("annotations.dir")
    )
@transform(
    "contigs.dir/*",
    regex("contigs.dir/(.*)\.fa"),
    r"annotations\.dir/\1.gff",
    r"\1"
    )
def prokka(infile, outfile, idd):
    
    ''' Annotate with prokka '''

    to_cluster = True

    statement = '''
    prokka --centre X --compliant %(infile)s --outdir annotations.dir --force --prefix %(idd)s
    '''

    P.run(statement)

# }}}
# roary {{{
@follows(
    prokka
    )
@transform(
    "annotations.dir",
    regex("annotations\.dir"),
    "roary.dir"
    )
def roary(infile, outfile):

    ''' Roary '''

    to_cluster = True

    statement = '''
    roary -f %(outfile)s -e -n -v -r %(infile)s/*.gff 
    '''

    P.run(statement)

    pass
# }}}
# distanceFromTree {{{
@follows(roary)
@transform(
    "roary.dir/*",
    suffix(".newick"),
    "distance.tsv"
    )
def distanceFromTree(infile, outfile):
    
    '''Get distances from a phylogeny tree that has been midpoint rooted'''

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
    mkdir("phenos.dir")
    )
@split(
    "phenos.tsv",
    "phenos.dir/*.tsv"
    )
def splitPhenos(infile, outfile):
    ''' split a tsv file into multiple tsv files by column '''

    statement = '''
    cols=`awk -F"\\t" '{print NF; exit}' %(infile)s` &&
    for col in $(seq 2 $cols); do
        pheno=`awk -F"\\t" -v col=$col 'NR==1{print $col}' %(infile)s` &&
        awk -F"\\t" -v col=$col '{print $1"\\t"$col}' %(infile)s > phenos.dir/${pheno}.tsv;
    done
    '''

    statement = '''
    echo %(infile)s &&
    echo %(outfile)s
    '''

    P.run(statement)

# }}}
# pyseer {{{
@follows(
    fsm,
    distanceFromTree,
    splitPhenos
    )
@transform(
    "phenos.dir/*",
    regex("phenos/(.*)\.tsv"),
    r"pyseer.dir/\1.assoc",
    extras = [distanceFromTree, fsm, r"\1"]
    )
def pyseer(infiles, outfile, idd):

    pheno = infiles[0][0]
    distance = infiles[1][0]
    kmers = infiles[1][1]

    statement = '''
    pyseer
        --phenotypes=%(pheno)s
        --kmers=%(kmers)s
        --distances=%(distance)s
        --min-af=0.01
        --max-af=0.99
        --cpu=15
        --filter-pvalue=1E-8
        > %(outfile)s
    '''

# }}}
# makeRefList {{{
@follows (
    prokka
    )
@transform(
    "prokka.dir",
    regex("prokka\.dir"),
    "ref.txt",
    add_inputs = ["ref.dir"]
    )
def makeRefList(infiles, outfile):

    ''' Make a list of references for kmer mapping '''

    prokka_dir = infiles[0][0]
    ref_dir = infiles[1][0]

    statement = '''
    ls %(prokka_dir)s | grep ".gff" | awk '{print $1"\t%(prokka_dir)/"$1"\tref" }' > %(outfile)s &&
    ls %(ref_dir)s  | awk '{print $1"\t%(ref_dir)/"$1"\tref" }' 
    '''
    P.run(statement)

# }}}
# mapKmers {{{
@follows(
    pyseer,
    makeRefList
    )
@transform(
    pyseer,
    regex("pyseer.dir/(.*).assoc"),
    r"map.dir/\1\.txt",
    add_inputs = [makeRefList]
    )
def mapKmers(infiles, outfile):

    infile = infiles[0][0]
    ref_genomes = infiles[1][0]

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    statement = '''
    python %(PY_SRC_PATH)/summarise_annotations.py
    '''
# }}}
# full {{{
@follows (
    mapKmers
    )
def full():
    pass

# }}}

def main():
    P.main(sys.argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
