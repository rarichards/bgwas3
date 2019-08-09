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
    [r"annotations/\1.gff", r"annotations/\1.tsv"],
    r"\1"
    )
def prokka(infile, outfile, idd):
    
    ''' Annotate with prokka '''

    gff = outfiles[0]
    tsv = outfiles[1]

    to_cluster = True

    statement = '''
    prokka --centre X --compliant %(infile)s --outdir extra/prokka --force --prefix %(idd)s &&
    mv extra/prokka/%(idd)s.gff %(gff)s
    mv extra/prokka/%(idd)s.tsv %(tsv)s
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
# qqplot {{{
@transform(
    pyseer,
    regex("^associations/(.*)_assoc.*$"),
    r"associations/\1_qqplot.png"
    )
def qqplot(infiles, outfile):

    assoc_gzip = infile[0][1]




# }}}
# makeRefList {{{
@follows(
    mkdir("refs")
    )
@merge(
    [prokka, [assembly], ["refs/*"]],
    "ref.txt"
    )
def makeRefList(infiles, outfile):

    ''' Make a list of references for kmer mapping '''

    temp = [str(i) for sublist in infiles for i in sublist]

    to_cluster = True

    gffs = [i for i in temp if re.match(r".*\.gff$", i)]
    refs = [i for i in gffs if re.match(r"^refs/.*", i)]
    drafts = [i for i in gffs if re.match(r"^annotations/.*", i)]

    with open(outfile, "w") as f:
        for gff in refs:
            idd = re.search("^.*/(.*)\.gff", gff).group(1)
            regex = r".*/" + idd + "\.(fa|fasta)"
            print(regex)
            fa = [i for i in temp if re.match(regex, i)][0]
            print(fa)
            f.write(fa + "\t" + gff + "\tref\n")
        for gff in drafts:
            idd = re.search("^.*/(.*)\.gff", gff).group(1)
            regex = ".*/" + idd + "\.(fa|fasta)"
            fa = [i for i in temp if re.match(regex, i)][0]
            f.write(fa + "\t" + gff + "\tdraft\n")

# }}}
# gff2tsv {{{
# @follows(
#     mkdir("refs")
#     )
# @transform(
#     "refs/*",
#     regex("refs/(.*)\.gff"),
#     r"refs/\1.tsv"
#     )
# def gff2tsv(infile, outfile):

#     R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

#     statement = '''
#     Rscript %(R_SRC_PATH)s/gff2tsv.R %(infile)s
#     '''
# # }}}
# bonferoni {{{
@transform(
    pyseer,
    regex("^associations/(.*)_assoc.*$"),
    [r"associations/\1_stats.txt", r"associations/\1_assoc_filtered.txt"]
    )
def bonferoni(infiles, outfiles):

    to_cluster = False

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    patterns = infiles[1]
    assoc_gzip = infiles[0]
    stats = outfiles[0]
    filtered = outfiles[1]

    statement = '''
    Rscript %(R_SRC_PATH)s/bonferoni.R %(patterns)s %(stats)s &&
    gzip -d -c %(assoc_gzip)s > temp.tsv &&
    wc -l temp.tsv | cut -f1 -d' ' | xargs -I @ echo -e 'kmers_tested\\t@' >> %(stats)s &&
    head -1 temp.tsv > %(filtered)s &&
    awk '$1=="bonf_thresh"{print $2}' %(stats)s | xargs -I @ sh -c 'awk '\\''$4<@{print $0}'\\'' temp.tsv' >> %(filtered)s &&
    wc -l %(filtered)s | cut -f1 -d' ' | xargs -I @ echo -e 'significant_kmers\\t@' >> %(stats)s &&
    rm temp.tsv
    '''

    P.run(statement)

# }}}
# seer2fa {{{
@transform(
    bonferoni,
    regex("^associations/(.*)_stats.*$"),
    r"associations/\1_assoc_filtered.fa"
    )
def seer2fa(infiles, outfile):

    to_cluster = False

    print(infiles)
    print(outfile)

    with open(outfile, "w") as file_out:
        with open(infiles[1]) as file_in:
            kmers_sum = 0
            header = file_in.readline()
            for kmer in file_in:
                kmers_sum +=1
                file_out.write(">" + str(kmers_sum) + "\n" + kmer.split("\t")[0] + "\n")

# }}}
# filter {{{
@transform(
    [pyseer, bonferoni],
    regex(r"^associations/(.*)_(assoc|stats).*$"),
    r"associations/\1_assoc_filtered.txt"
    )
def filter(infiles, outfile):

    statement = '''
    cat <(head -1 penicillin_kmers.txt) <(awk '$4<1.90E-08 {print $0}' penicillin_kmers.txt) > significant_kmers.txt
    '''
    print(infiles)
    print(outfile)
# }}}
# bwa {{{
@transform(
    assembly,
    suffix(".fa"),
    [".amb", ".ann", ".bwt", ".pac", ".sa"]
    )
def bwa(infile, outfiles):
    
    statement = '''
    bwa index %(infile)s
    '''

    P.run(statement)

# }}}
# mapKmers {{{
@transform(
    bonferoni,
    regex(r"^associations/(.*)_stats.txt$"),
    add_inputs(makeRefList),
    r"associations/\1_map.txt"
    )
def mapKmers(infiles, outfile):

    to_cluster = False

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))
    
    filtered = infiles[0][1]
    ref_list = infiles[1]

    statement = '''
    python %(PY_SRC_PATH)s/annotate_kmers.py %(filtered)s %(ref_list)s %(outfile)s
    '''

    P.run(statement)

# }}}
# countGeneHits {{{
@transform(
    mapKmers,
    regex("^associations/(.*)_map.txt$"),
    r"associations/\1_hits.txt"
    )
def countGeneHits(infile, outfile):

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    statement = '''
    python %(PY_SRC_PATH)s/summarise_annotations.py %(infile)s > %(outfile)s
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
