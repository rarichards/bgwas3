import os
import sys
from ruffus import *
import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as iotools

import re

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
    to_cluster = False
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
    fsm-lite -l ../contigs_list.txt -s 6 -S 610 -v -t fsm_kmers | gzip -c > ../%(outfile)s
    '''

    P.run(statement)

# }}}
# prokka {{{
@follows(
    mkdir("extra/prokka"),
    mkdir("annotations")
    )
@transform(
    assembly,
    regex("contigs/(.*)\.fa"),
    r"annotations/\1.gff",
    r"\1"
    )
def prokka(infile, outfile, idd):
    
    ''' Annotate with prokka '''

    to_cluster = True

    statement = '''
    prokka --centre X --compliant %(infile)s --outdir extra/prokka --force --prefix %(idd)s &&
    mv extra/prokka/%(idd)s.gff %(outfile)s
    '''

    P.run(statement)

# }}}
# roary {{{
@follows(
    prokka,
    mkdir("extra/roary")
    )
@merge(
    prokka,
    "tree.newick"
    )
def roary(infile, outfile):

    ''' Make tree with Roary '''

    to_cluster = True

    statement = '''
    roary -f extra/roary -e -n -v -r annotations/*.gff &&
    cp extra/roary/accessory_binary_genes.fa.newick %(outfile)s
    '''

    P.run(statement)

    pass

# }}}
# distanceFromTree {{{
@transform(
    roary,
    regex("roary/accessory_binary_genes\.fa\.newick"),
    "distances.tsv"
    )
def distanceFromTree(infile, outfile):
    
    ''' Get distances from a phylogeny tree that has been midpoint rooted '''

    to_cluster = False

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

    ''' Split a tsv file into multiple tsv files by column '''

    to_cluster = False

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
    mkdir("associations"),
    mkdir("extra/pyseer")
    )
@transform(
    splitPhenos,
    regex("phenos/(.*)\.tsv"),
    add_inputs(distanceFromTree, fsm),
    r"associations/\1.assoc.gz",
    r"\1"
    )
def pyseer(infiles, outfile, idd):

    to_cluster = True

    pheno = infiles[0]
    distances = infiles[1]
    kmers = infiles[2]

    print(pheno)
    print(distances)
    print(kmers)

    statement = '''
    pyseer 
        --lmm 
        --phenotypes %(pheno)s
        --kmers %(kmers)s
        --similarity %(distances)s
        --output-patterns extra/pyseer/%(idd)s_patterns.txt
        --cpu 8 
        | gzip -c 
        > %(outfile)s
    '''

    P.run(statement)

# }}}
# makeRefList {{{
@follows(
    mkdir("refs")
    )
@merge(
    [prokka, assembly, "refs/*"],
    "ref.txt"
    )
def makeRefList(infiles, outfile):

    ''' Make a list of references for kmer mapping '''

    to_cluster = True

    gffs = list(filter(re.compile(".*\.gff$").match, infiles))
    refs = list(filter(re.compile("refs/.*").match, gffs))
    drafts = list(filter(re.compile("annotations/.*").match, gffs))

    with open(outfile, "w") as f:
        for gff in refs:
            idd = re.search("^.*/(.*)\.gff", gff).group(1)
            regex = ".*/" + idd + "\.(fa|fasta)"
            fa = list(filter(re.compile(regex).match, infiles))[0]
            f.write(fa + "\t" + gff + "\tref\n")
        for gff in drafts:
            idd = re.search("^.*/(.*)\.gff", gff).group(1)
            regex = ".*/" + idd + "\.(fa|fasta)"
            fa = list(filter(re.compile(regex).match, infiles))[0]
            f.write(fa + "\t" + gff + "\tdraft\n")


# }}}
# mapKmers {{{
@follows(
    mkdir("maps")
    )
@transform(
    pyseer,
    regex(r"^associations/(.*)\.assoc$"),
    add_inputs(makeRefList),
    r"maps/\1.txt"
    )
def mapKmers(infiles, outfile):

    to_cluster = True

    kmers = infiles[0]
    refs = infiles[1]

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    statement = '''
    python %(PY_SRC_PATH)s/annotate_kmers.py %(kmers)s %(refs)s %(outfile)s
    '''

    P.run(statement)

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
