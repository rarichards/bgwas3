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

# makeDirs {{{
@originate(
    "results"
    )
def makeDirs(outfile):

    ''' Copy plot directory (index.html) to working directory '''

    output_dir = os.path.dirname(os.path.realpath(__file__)) + "/plot"
    cwd = os.getcwd()
    shutil.copytree(output_dir, os.getcwd())

# }}}
# assembly {{{
#@transform(
#    "fastqs/*",
#    regex(".*"),
#    "contigs"
#    )
@originate(
    "contigs"
    )
def assembly(infile, outfile):

    ''' Contig assembly '''

    # TODO this step with options from yml

    to_cluster = False
    pass

# }}}
# mineKmers {{{
@merge(
    assembly,
    "kmers.txt.gz"
    )
def mineKmers(infile, outfile):

    ''' Kmer mining/ counting with fsm-lite (make a gzipped kmers file) '''

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

    ''' Annotate genomes with prokka (generate gff file for each contig) '''

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

    ''' Perform pangenome analalysis with Roary (based on difference in gene 
    content) and generate a phylogenetic tree (.newick) file for use in mixed 
    effects association testing '''

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
    
    ''' Get distances from a phylogenetic tree that has been midpoint 
    rooted '''

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

    ''' Plot phylogenetic trees with filters based on clade and phenotype '''

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

    ''' Split the main tsv file phenotype columns into their own tsv files 
    (to be used in seperated association tests) '''

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
    mkdir("results"),
    )
@transform(
    splitPhenos,
    regex("phenos/(.*)\.tsv"),
    add_inputs(distanceFromPangenome, mineKmers),
    [r"results/\1_assocs.txt.gz", r"results/\1_patterns.txt"],
    r"\1"
    )
def testAssociations(infiles, outfiles, idd):

    ''' Association test kmers and phenotypes with pyseer '''

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
    regex("^results/(.*)_assocs.*$"),
    [r"results/plot/p/\1_hist.png", r"results/plot/p/\1_qq.png"],
    r"\1",
    "output/p"
    )
def qqplot(infiles, outfile, pheno, outdir):

    ''' Plot p-values (histogram and qqlot) for each association test '''

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
    regex("^results/(.*)_assocs.txt.gz"),
    [r"results/\1_stats.txt", r"results/\1_assocs_filtered.txt"]
    )
def bonferoniFilter(infiles, outfiles):

    ''' Filter kmers from the output of association based on their p-value,
    using bonferoni method to set a threshold '''

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
    wc -l %(filtered)s | awk '{sum = $1 - 1; print sum}' | xargs -I @ echo -e 'significant_kmers\\t@' >> %(stats)s &&
    rm temp.tsv
    '''

    P.run(statement)

# }}}

# annotation2bed {{{
@transform(
    annotateGenomes,
    regex("annotations/(.*)\.gff"),
    r"annotations/\1.bed"
    )
def annotation2bed(infile, outfile):

    ''' Convert gff files into bed files. Also remove annotations that don't 
    correspind to named genes (eg. hypothetical proteins) '''

    if os.stat(infile).st_size == 0:
        statement = '''
        touch %(outfile)s
        '''
    else:
        statement = '''
        gff2bed < %(infile)s | grep "gene=" > %(outfile)s
        '''
    P.run(statement, to_cluster=False)

# }}}
# ref2bed {{{
@transform(
    "refs/*",
    regex("refs/(.*)\.gff*"),
    r"refs/\1.bed"
    )
def ref2bed(infile, outfile):

    ''' Convert reference gff files into bed files. Also merge overlapping
    annotations (eg 'CDS' and 'Gene' entries that are the same, but have different 
    information. Also remove 'Region' entries. '''

    statement = '''
    gff2bed < %(infile)s | awk '$8 != "region"{print $0}' | bedtools merge -c 10 -o collapse -delim "|" > %(outfile)s
    '''

    P.run(statement, to_cluster=False)

# }}}
# makeRefList {{{
@merge(
    [[ref2bed],[annotation2bed],[assembly, "refs/*.fa"]],
    "ref.txt"
    )
def makeRefList(infiles, outfile):

    ''' Make a list of references for kmer mapping (references first, then 
    draft annotations) '''

    print(infiles)

    ref_beds = infiles[0]
    annotation_beds = infiles[1]
    fas = infiles[2]

    with open(outfile, "w") as f:
        for bed in ref_beds:
            idd = re.search("^.*/(.*)\.bed", bed).group(1)
            regex = r".*/" + idd + "\.fa"
            fa = [i for i in fas if re.match(regex, i)][0]
            f.write(fa + "\t" + bed + "\n")
        for bed in annotation_beds:
            if os.stat(bed).st_size != 0:
                idd = re.search("^.*/(.*)\.bed", bed).group(1)
                regex = r".*/" + idd + "\.fa"
                fa = [i for i in fas if re.match(regex, i)][0]
                f.write(fa + "\t" + bed + "\n")


# }}}
# bwaIndex {{{
@transform(
    [assembly, "refs/*"],
    suffix(".fa"),
    [".fa.amb", ".fa.ann", ".fa.bwt", ".fa.pac", ".fa.sa"]
    )
def bwaIndex(infile, outfiles):
    
    ''' BWA index (for Kmer mapping) '''
    
    statement = '''
    bwa index %(infile)s
    '''

    P.run(statement)

# }}}
# mapKmers{{{
@transform(
    bonferoniFilter,
    regex(r"^(.*)_stats.txt$"),
    add_inputs(makeRefList, ref2bed, annotation2bed, bwaIndex),
    [r"\1_maps.txt", r"\1_gene_info.txt"],
    r"\1"
    )
def mapKmers(infiles, outfiles, prefix):

    ''' Map kmers to annotation files '''

    to_cluster = False

    PY_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "python"))

    filtered = infiles[0][1]
    ref_list = infiles[1]
    print(prefix + filtered)

    statement = '''
    python %(PY_SRC_PATH)s/mapKmers.py %(filtered)s %(ref_list)s --prefix %(prefix)s
    '''

    P.run(statement, to_cluster=False)

# }}}

# summariseGenes {{{
@transform(
    mapKmers,
    regex("^results/(.*)_map.txt$"),
    r"results/\1_genes.txt"
    )
def countGeneHits(infiles, outfile):

    ''' Count the number of kmer hits per gene, and calculate other
    statistics '''

    maps = infiles[0]
    gene_info = infiles[1]

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
    
    # TODO

    pass

# }}}

# plotGenes {{{
@merge(
    [countGeneHits, pathwayAnalysis],
    "results/plot/dat.js"
    )
def visualise(infile, outfile):

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
