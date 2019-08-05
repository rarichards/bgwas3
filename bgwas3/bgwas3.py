import os
import sys
import shutil
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
    "kmers.txt.gz"
    )
def fsm(infile, outfile):

    ''' Kmer mining/ counting with fsm-lite '''

    print(PARAMS)

    to_cluster = True

    statement = '''
    ls contigs | awk -F. '{print $1 "\t" $0}' > contigs_list.txt &&
    cd contigs &&
    fsm-lite 
        -l ../contigs_list.txt 
        -m %(fsm_kmer-min)s
        -M %(fsm_kmer-max)s
        -v
        -t kmers
        | gzip -c > ../%(outfile)s
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
    regex("tree\.newick"),
    "distances.tsv"
    )
def distanceFromTree(infile, outfile):
    
    ''' Get distances from a phylogeny tree that has been midpoint rooted '''

    to_cluster = False

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    print("\n" + infile + "\n")
    print("\n" + outfile + "\n")

    statement = '''
    python %(PY_SRC_PATH)s/phylogeny_distance.py 
        --calc-C %(infile)s 
        > %(outfile)s
    '''
    
    P.run(statement)

# }}}
# plotTrees {{{
@follows(
    mkdir("out"),
    mkdir("out/static"),
    mkdir("out/static/trees")
    )
@split(
    [roary, "phenos.tsv"],
    "out/static/trees/*.png"
    )
def plotTrees(infiles, outfiles):

    tree = infiles[0]
    phenos = infiles[1]

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement =  '''
    Rscript %(R_SRC_PATH)s/plotTrees.R trees.newick phenos.tsv out/static/trees
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

    ''' Split the main tsv file phenotype columns into their own tsv files '''

    to_cluster = False

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement = '''
    Rscript %(R_SRC_PATH)s/splitPhenos.R %(infile)s phenos
    '''

    # statement = '''
    # cols=`awk -F"\\t" '{print NF; exit}' %(infile)s` &&
    # for col in $(seq 2 $cols); do
    #     pheno=`awk -F"\\t" -v col=$col 'NR==1{print $col}' %(infile)s` &&
    #     awk -F"\\t" -v col=$col '{print $1"\\t"$col}' %(infile)s > phenos/${pheno}.tsv;
    # done
    # '''

    P.run(statement)

# }}}
# pyseer {{{
@follows(
    mkdir("associations"),
    )
@transform(
    splitPhenos,
    regex("phenos/(.*)\.tsv"),
    add_inputs(distanceFromTree, fsm),
    [r"associations/\1_assoc.txt.gz", r"associations/\1_patterns.txt"],
    r"\1"
    )
def pyseer(infiles, outfiles, idd):

    to_cluster = True

    pheno = infiles[0]
    distances = infiles[1]
    kmers = infiles[2]

    assoc = outfiles[0]
    patterns = outfiles[1]

    print(pheno)
    print(distances)
    print(kmers)

    statement = '''
    pyseer 
        --lmm 
        --phenotypes %(pheno)s
        --kmers %(kmers)s
        --similarity %(distances)s
        --output-patterns %(patterns)s
        --cpu 8 
        | gzip -c 
        > %(assoc)s
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
# gff2tsv {{{
@follows(
    mkdir("refs")
    )
@transform(
    "refs/*",
    regex("refs/(.*)\.gff"),
    r"refs/\1.tsv"
    )
def gff2tsv(infile, outfile):

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement = '''
    Rscript %(R_SRC_PATH)s/gff2tsv.R %(infile)s
    '''
# }}}
# bonferoni {{{
@transform(
    pyseer,
    regex(r"^associations/(.*)_patterns.txt$"),
    r"assocations/\1_stats.txt"
    )
def bonferoni(infiles, outfile):

    to_cluster = False

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    print(infiles)

    assoc_gzip = infiles[0]
    assoc = infiles[0][:-3]

    statement = '''
    gzip -d %(assoc_gzip)s &&
    Rscript %(R_SRC_PATH)s/bonferoni.R %(assoc)s > %(outfile)s
    gzip %(assoc)s
    '''

    P.run(statement)

# }}}
# filter {{{
@transform(
    [pyseer, bonferoni],
    regex(r"^associations/(.*)_(assoc|stats).*$"),
    r"associations/\1_assoc_filtered.txt"
    )
def filter(infiles, outfile):
    print(infiles)
    print(outfiles)
# }}}
# mapKmers {{{
@follows(
    mkdir("maps")
    )
@transform(
    filter,
    regex(r"^associations/(.*)_assoc\.txt\.gz$"),
    add_inputs(makeRefList),
    r"maps/\1_map.txt.gz"
    )
def mapKmers(infiles, outfile):

    to_cluster = True

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    assoc_gzip = infiles[0]
    assoc = infiles[0][:-3]
    ref_list = infiles[1]
    maps = outfile[:-3]
    
    reflist = infiles[1]

    statement = '''
    gzip -d %(assoc_gzip)s &&
    python %(PY_SRC_PATH)s/annotate_kmers.py %(assoc)s %(ref_list)s %(maps)s &&
    gzip %(assoc)s &&
    gzip %(maps)s
    '''

    print(statement)

    P.run(statement)

# }}}
# countGeneHits {{{
@transform(
    filter,
    regex("^maps/(.*)_maps_filtered.txt$"),
    r"maps/\1_hits.txt"
    )
def countGeneHits(infile, outfile):

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    infile_unzip = infile[:-3]
    outfile_unzip = outfile[:-3]

    statement = '''
    gzip -d %(infile)s &&
    python %(PY_SRC_PATH)s/summarise_annotations.py %(infile_unzip)s > %(outfile_unzip)s &&
    gzip %(outfile_unzip)s
    '''

    P.run(statement)

# }}}
# pathwayAnalysis {{{
@transform(
    countGeneHits,
    regex("^maps/(.*)_hits.txt$"),
    r"maps/\1_pathways.tsv"
    )
def pathwayAnalysis(infiles, outfile):
    pass

# }}}
# visualise {{{
@follows(
    mkdir("viz")
    )
@merge(
    [countGeneHits, pathwayAnalysis],
    "visual"
    )
def visualise(infile, outfile):
    pass

# }}}
# full {{{
@follows (
    visualise
    )
def full():
    pass

# }}}

def main():
    P.main(sys.argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
