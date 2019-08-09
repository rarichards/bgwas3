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

# generateOutputDir {{{
@originate(
    "output"
    )
def generateOutputDir(outfile):
    output_dir = os.path.dirname(os.path.realpath(__file__)) + "/output"
    cwd = os.getcwd()
    shutil.copytree(output_dir, os.getcwd())

# }}}
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
# mineKmers {{{
@merge(
    assembly,
    "kmers.txt.gz"
    )
def mineKmers(infile, outfile):

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
# annotateGenomes {{{
@follows(
    mkdir("annotations")
    )
@transform(
    assembly,
    regex("contigs/(.*)\.fa"),
    r"annotations/\1.gff",
    r"\1"
    )
def annotateGenomes(infile, outfile, idd):
    ''' Annotate genomes with prokka '''

    statement = '''
    prokka --centre X --compliant %(infile)s --outdir annotations --force --prefix %(idd)s
    '''

    P.run(statement, to_cluster=True)

# }}}
# pangenomeAnalysis {{{
@merge(
    annotateGenomes,
    "pangenome/accessory_binary_genes.fa.newick"
    )
def pangenomeAnalysis(infile, outfile):

    ''' Make tree with Roary '''

    os.mkdir("pangenome")
    statement = '''
    roary -f pangenome -e -n -v -r annotations/*.gff
    '''

    P.run(statement, to_cluster=True)

# }}}
# distancesFromPangenome {{{
@transform(
    pangenomeAnalysis,
    regex("pangenome/accessory_binary_genes.fa.newick"),
    "distances.tsv"
    )
def distanceFromPangenome(infile, outfile):
    
    ''' Get distances from a phylogeny tree that has been midpoint rooted '''

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    newick = infile + "/accessory_binary_genes.fa.newick"

    statement = '''
    python %(PY_SRC_PATH)s/phylogeny_distance.py 
        --calc-C %(newick)s 
        > %(outfile)s
    '''
    
    P.run(statement, to_cluster=False)

# }}}
# plotPhyloTree {{{
@follows(
    generateOutputDir
    )
@split(
    [pangenomeAnalysis, "phenos.tsv"],
    "output/trees/*"
    )
def plotTrees(infiles, outfiles):

    tree = "pangenome/accessory_binary_genes.fa.newick"
    phenos = infiles[1]

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement =  '''
    Rscript %(R_SRC_PATH)s/plotTrees.R trees.newick phenos.tsv output/trees
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
# testAssociations {{{
@follows(
    mkdir("associations"),
    )
@transform(
    splitPhenos,
    regex("phenos/(.*)\.tsv"),
    add_inputs(distanceFromPangenome, mineKmers),
    [r"associations/\1_assoc.txt.gz", r"associations/\1_patterns.txt"],
    r"\1"
    )
def testAssociations(infiles, outfiles, idd):

    pheno = infiles[0]
    distances = infiles[1]
    kmers = infiles[2]

    assoc = outfiles[0]
    patterns = outfiles[1]

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

    P.run(statement, to_cluster=True)

# }}}
# plotP {{{
@follows(
    generateOutputDir
    )
@transform(
    testAssociations,
    regex("^associations/(.*)_assoc.*$"),
    [r"output/p/\1_hist.png", r"output/p/\1_qq.png"],
    r"\1",
    "output/p"
    )
def qqplot(infiles, outfile, pheno, outdir):

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement =  '''
    zcat %(infile)s | awk '{print $4}' > temp &&
    Rscript %(R_SRC_PATH)s/plotP.R temp --out %(outdir)s --prefix %(pheno)s &&
    rm temp
    '''

    P.run(statement)

# }}}
# bonferoniFilter {{{
@transform(
    testAssociations,
    regex("^associations/(.*)_assoc.*$"),
    [r"associations/\1_stats.txt", r"associations/\1_assoc_filtered.txt"]
    )
def bonferoniFilter(infiles, outfiles):

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
# filterAnnotations {{{
@transform(
    annotateGenomes,
    regex("annotations/(.*)\.gff"),
    r"annotations/\1_genes_only.gff"
    )
def filterAnnotations(infile, outfile):
    
    statement = '''
    cat %(infile)s | grep "gene=" > %(outfile)s
    '''

    P.run(statement, to_cluster=False)

# }}}
# makeRefList {{{
@merge(
    [[filterAnnotations], [assembly]],
    "ref.txt"
    )
def makeRefList(infiles, outfile):

    ''' Make a list of references for kmer mapping '''

    gffs = infiles[0]
    fas = infiles[1]

    with open(outfile, "w") as f:
        for gff in gffs:
            idd = re.search("^.*/(.*)\_genes_only.gff", gff).group(1)
            regex = r".*/" + idd + "\.(fa|fasta)"
            fa = [i for i in fas if re.match(regex, i)][0]
            f.write(fa + "\t" + gff + "\tdraft\n")

# }}}
# mapKmers{{{
@transform(
    bonferoniFilter,
    regex(r"^associations/(.*)_stats.txt$"),
    add_inputs(makeRefList),
    r"associations/\1_maps.txt"
    )
def mapKmers(infiles, outfile):

    to_cluster = False

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))
    
    filtered = infiles[0][1]
    ref_list = infiles[1]

    statement = '''
    python %(PY_SRC_PATH)s/annotate_kmers.py %(filtered)s %(ref_list)s %(outfile)s
    '''

    P.run(statement, to_cluster=False)

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
