import os
import sys
import shutil
import re
from ruffus import *
from cgatcore import pipeline as P

PARAMS = P.get_parameters([
    "%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "../pipeline.yml",
    "pipeline.yml"
])

# make_dirs {{{
@follows(
    mkdir("results")
)
@originate(
    "results/plots"
)
def make_dirs(outfile):

    ''' Copy plot directory (index.html) to working directory '''

    output_dir = os.path.dirname(os.path.realpath(__file__)) + "/plots"
    shutil.copytree(output_dir, os.getcwd() + "/" + outfile)

# }}}
# assembly {{{
# @transform(
#    "fastqs/*",
#    regex(".*"),
#    "contigs"
#    )

@originate(
    "contigs"
)
def assembly(outfile):

    ''' Contig assembly '''

    # TODO this step with options from yml

# }}}
# mine_kmers {{{
@merge(
    assembly,
    "kmers.txt.gz"
)
def mine_kmers(infile, outfile):

    ''' Kmer mining/ counting with fsm-lite (make a gzipped kmers file) '''

    print(PARAMS)

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

    P.run(statement, to_cluster=True)

# }}}
# annotate {{{
@follows(
    mkdir("annotations")
)
@transform(
    assembly,
    regex("contigs/(.*).fa"),
    r"annotations/\1.gff",
    r"\1"
)
def annotate(infile, outfile, idd):

    ''' Annotate genomes with prokka (generate gff file for each contig) '''

    statement = '''
    prokka --centre X --compliant %(infile)s --outdir annotations --force
                                                         --prefix %(idd)s
    '''

    P.run(statement, to_cluster=True)

# }}}
# pangenome_analysis {{{
@merge(
    annotate,
    "pangenome/accessory_binary_genes.fa.newick"
)
def pangenome_analysis(infile, outfile):

    ''' Perform pangenome analalysis with Roary (based on difference in gene
    content) and generate a phylogenetic tree (.newick) file for use in mixed
     effects association testing '''

    os.mkdir("pangenome")
    statement = '''
    roary -f pangenome -e -n -v -r annotations/*.gff
    '''

    P.run(statement, to_cluster=True)

# }}}
# distance_from_tree {{{
@transform(
    pangenome_analysis,
    regex("pangenome/accessory_binary_genes.fa.newick"),
    "distances.tsv"
)
def distance_from_tree(infile, outfile):

    ''' Get distances from a phylogenetic tree that has been midpoint
    rooted '''

    PY_SRC_PATH = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "python")
    )

    newick = infile + "/accessory_binary_genes.fa.newick"

    statement = '''
    python %(PY_SRC_PATH)s/distance_from_tree.py
        --calc-C %(newick)s
        > %(outfile)s
    '''

    P.run(statement, to_cluster=False)

# }}}

# plot_trees {{{
@follows(
    make_dirs
)
@split(
    [pangenome_analysis, "phenos.tsv"],
    "output/trees/*"
)
def plot_trees(infiles, outfiles):

    ''' Plot phylogenetic trees with filters based on clade and phenotype '''

    tree = "pangenome/accessory_binary_genes.fa.newick"
    phenos = infiles[1]

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement = '''
    Rscript %(R_SRC_PATH)s/plot_trees.R %(tree)s %(phenos)s results/plots
    '''

    P.run(statement)

# }}}
# split_phenos {{{
@follows(
    mkdir("phenos")
)
@split(
    "phenos.tsv",
    "phenos/*.tsv"
)
def split_phenos(infile, outfiles):

    ''' Split the main tsv file phenotype columns into their own tsv files
    (to be used in seperated association tests) '''


    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement = '''
    Rscript %(R_SRC_PATH)s/split_phenos.R %(infile)s phenos
    '''

    P.run(statement, to_cluster=False)

# }}}

# test_assoc {{{
@follows(
    mkdir("results"),
)
@transform(
    split_phenos,
    regex("phenos/(.*).tsv"),
    add_inputs(distance_from_tree, mine_kmers),
    [r"results/\1_assocs.txt.gz", r"results/\1_patterns.txt"],
    r"\1"
)
def test_assoc(infiles, outfiles, idd):

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
# plot_ps {{{
@follows(
    make_dirs
)
@transform(
    test_assoc,
    regex("^results/(.*)_assocs.*$"),
    [r"results/plot/p/\1_hist.png", r"results/plot/p/\1_qq.png"],
    r"\1",
    "output/p"
)
def plot_ps(infiles, outfile, pheno, outdir):

    ''' Plot p-values (histogram and qqplot) for each association test '''

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement = '''
    zcat %(infile)s | awk '{print $4}' > temp &&
    Rscript %(R_SRC_PATH)s/plot_ps.R temp --out %(outdir)s --prefix %(pheno)s &&
    rm temp
    '''

    P.run(statement)

# }}}
# bonferoni{{{
@transform(
    test_assoc,
    regex("^results/(.*)_assocs.txt.gz"),
    [r"results/\1_stats.txt", r"results/\1_assocs_filtered.txt"]
)
def bonferoni(infiles, outfiles):

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
    wc -l temp.tsv | cut -f1 -d' ' | xargs -I @ echo -e 'kmers_tested\\t@'
    >> %(stats)s &&
    head -1 temp.tsv > %(filtered)s &&
    awk '$1=="bonf_thresh"{print $2}' %(stats)s |
    xargs -I @ sh -c 'awk '\\''$4<@{print $0}'\\'' temp.tsv' >> %(filtered)s &&
    wc -l %(filtered)s | awk '{sum = $1 - 1; print sum}'
    | xargs -I @ echo -e 'significant_kmers\\t@' >> %(stats)s &&
    rm temp.tsv
    '''

    P.run(statement)

# }}}

# annotation2bed {{{
@transform(
    annotate,
    regex("annotations/(.*).gff"),
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
    regex("refs/(.*).gff*"),
    r"refs/\1.bed"
)
def ref2bed(infile, outfile):

    ''' Convert reference gff files into bed files. Also merge overlapping
    annotations (eg 'CDS' and 'Gene' entries that are the same, but have
    different information. Also remove 'Region' entries. '''

    statement = '''
    gff2bed < %(infile)s | awk '$8 != "region"{print $0}' |
    bedtools merge -c 10 -o collapse -delim "|" > %(outfile)s
    '''

    P.run(statement, to_cluster=False)

# }}}
# make_ref_list {{{
@merge(
    [[ref2bed], [annotation2bed], [assembly, "refs/*.fa"]],
    "ref.txt"
)
def make_ref_list(infiles, outfile):

    ''' Make a list of references for kmer mapping (references first, then
    draft annotations) '''

    print(infiles)

    ref_beds = infiles[0]
    annotation_beds = infiles[1]
    fas = infiles[2]

    with open(outfile, "w") as f:
        for bed in ref_beds:
            idd = re.search("^.*/(.*).bed", bed).group(1)
            regex = r".*/" + idd + ".fa"
            fa = [i for i in fas if re.match(regex, i)][0]
            f.write(fa + "\t" + bed + "\n")
        for bed in annotation_beds:
            if os.stat(bed).st_size != 0:
                idd = re.search("^.*/(.*).bed", bed).group(1)
                regex = r".*/" + idd + ".fa"
                fa = [i for i in fas if re.match(regex, i)][0]
                f.write(fa + "\t" + bed + "\n")


# }}}
# bwa_index {{{
@transform(
    [assembly, "refs/*"],
    suffix(".fa"),
    [".fa.amb", ".fa.ann", ".fa.bwt", ".fa.pac", ".fa.sa"]
)
def bwa_index(infile, outfiles):

    ''' BWA index (for Kmer mapping) '''

    statement = '''
    bwa index %(infile)s
    '''

    P.run(statement)

# }}}
# map_kmers{{{
@transform(
    bonferoni,
    regex(r"^(.*)_stats.txt$"),
    add_inputs(make_ref_list, ref2bed, annotation2bed, bwa_index),
    [r"\1_maps.txt", r"\1_gene_info.txt"],
    r"\1"
)
def map_kmers(infiles, outfiles, prefix):

    ''' Map kmers to annotation files '''

    to_cluster = False

    PY_SRC_PATH = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "python")
    )

    filtered = infiles[0][1]
    ref_list = infiles[1]
    print(prefix + filtered)

    statement = '''
    python %(PY_SRC_PATH)s/map_kmers.py %(filtered)s %(ref_list)s
        --prefix %(prefix)s
    '''

    P.run(statement, to_cluster=False)

# }}}

# summarise_genes {{{
@transform(
    map_kmers,
    regex("^results/(.*)_map.txt$"),
    r"results/\1_genes.txt"
)
def summarise_genes(infiles, outfile):

    ''' Count the number of kmer hits per gene, and calculate other
    statistics '''

    maps = infiles[0]
    gene_info = infiles[1]

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement = '''
    python %(R_SRC_PATH)s/summarise_genes.R %(infile)s > %(outfile)s
    '''

    P.run(statement)

# }}}
# pathway_analysis {{{
@transform(
    summarise_genes,
    regex("^maps/(.*)_hits.txt$"),
    r"maps/\1_pathways.tsv"
)
def pathwayAnalysis(infiles, outfile):
    pass

# }}}

# plot_genes {{{
@merge(
    [summarise_genes],
    "results/plot/dat.js"
)
def visualise(infile, outfile):
    pass

# }}}

# full {{{

@follows(
    visualise
)
def full():
    pass

# }}}

def main():
    P.main(sys.argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
