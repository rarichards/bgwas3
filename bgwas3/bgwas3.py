import os
import sys
import shutil
import re
from ruffus import * 
from cgatcore import pipeline as P
from mako.template import Template
import json

'''
'''

PARAMS = P.get_parameters([
    "%s/pipeline.yml" % os.path.splitext(__file__)[0],
    "../pipeline.yml",
    "pipeline.yml"
])

# assembly {{{
@follows(
    mkdir("fastqs")
    )
@split(
    "fastqs",
    "contigs/*.fa"
)
def assembly(outfile):

    ''' 
    Contig assembly
    '''

# }}}
# mine_kmers {{{
@merge(
    assembly,
    "kmers.txt.gz"
)
def mine_kmers(infile, outfile):

    ''' Kmer mining and counting with fsm-lite
    
    :param infile: directory of geneomes (fasta files)
    :param outfile: gzipped file of Kmer patterns and genomes they are found in

    '''

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

# mash {{{
@merge(
    assembly,
    regex("^contigs/(.*).fa$"),
    "mash.tsv"
    )
def mash(infiles, outfile):

    '''
    '''

    statement = '''
    mash sketch -s 10000 -o mash_sketch assemblies/*.fa
    '''

    P.run(statement, to_cluster = False)

# }}}

# plot_trees {{{
@split(
    [pangenome_analysis, "phenos.tsv"],
    "results/plots/tree.html"
)
def plot_trees(infiles, outfile):

    ''' Plot phylogenetic trees with filters based on clade and phenotype '''

    template_path = os.path.dirname(os.path.realpath(__file__)) + "/plots/tree.html"
    template = Template(filename=template_path)

    json_phenos = ""
    with open(infiles[1], "r") as dat:
        titles = dat.readline().split(sep="\t")
        for l in dat:
            d = {}
            for t, f in zip(titles, l.split(sep="\t")):
                d[t.rstrip()] = f.rstrip()
            json_phenos += json.dumps(d, indent=4) + ","

    with open(infiles[0], "r") as tree_file:
        newick = tree_file.readline();

    page = template.render(json_phenos=json_phenos, newick=newick)
    
    with open(outfile, "w") as html_file:
        html_file.write(page)
   
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
# plot_phenos {{{
@follows(
    mkdir("results/plots"),
)
@transform(
    split_phenos,
    regex("^phenos/(.*).tsv$"),
    r"results/plots/\1_density.png",
    r"\1"
)
def plot_phenos(infile, outfile, pheno):

    ''' Density plot phenotypes '''

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

        
    statement = '''
    Rscript %(R_SRC_PATH)s/plot_density.R %(infile)s %(outfile)s --width 1 --height 1 --axis 0 --column 2
    '''

    P.run(statement)

# }}}

# test_assoc {{{
@follows(
    mkdir("results"),
)
@transform(
    split_phenos,
    regex("phenos/(.*).tsv"),
    add_inputs(distance_from_tree, mine_kmers),
    [r"results/\1_assocs.tsv.gz", r"results/patterns.txt"],
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
    mkdir("results/plots"),
    mkdir("results/plots/p")
)
@transform(
    test_assoc,
    regex("^results/(.*)_assocs.tsv.gz$"),
    [r"results/plots/p/\1_hist.png", r"results/plots/p/\1_qq.png", r"results/plots/p/\1_unadj_hist.png", r"results/plots/p/\1_unadj_qq.png", r"results/plots/\1_p.html"],
    r"\1"
)
def plot_ps(infiles, outfiles, pheno):

    ''' Plot p-values (histogram and qqplot) for each association test '''

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    PY_SRC_PATH = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "python")
    )

    assoc = infiles[0]

    p_hist = outfiles[0]
    p_qq = outfiles[1]
    p_unadj_hist = outfiles[2]
    p_unadj_qq = outfiles[3]

    template_path = os.path.dirname(os.path.realpath(__file__)) + "/plots/p.html"
    template = Template(filename=template_path)

    html_file = open(outfiles[4])
    page = template.render(pheno=pheno, p_hist=p_hist, p_qq=p_qq)
    html_file.write(page)
    html_file.close()
    
    statement = '''
    zcat %(assoc)s | awk '{print $4}' > %(pheno)s_temp_p &&
    zcat %(assoc)s | awk '{print $3}' > %(pheno)s_temp_p_unadj &&
    Rscript %(R_SRC_PATH)s/plot_ps.R %(pheno)s_temp_p --output %(p_hist)s &&
    Rscript %(R_SRC_PATH)s/plot_ps.R %(pheno)s_temp_p_unadj --output %(p_unadj_hist)s &&
    python %(PY_SRC_PATH)s/plot_qq.py %(pheno)s_temp_p --output %(p_qq)s &&
    python %(PY_SRC_PATH)s/plot_qq.py %(pheno)s_temp_p_unadj --output %(p_unadj_qq)s_anadj &&
    rm %(pheno)s_temp_p*
    '''

    P.run(statement)

# }}}

# bonferoni {{{
@merge(
    test_assoc,
    "results/bonf.txt"
)
def bonferoni(infiles, outfile):

    ''' Filter kmers from the output of association based on their p-value,
    using bonferoni method to set a threshold '''

    to_cluster = False

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    patterns = infiles[0][1]

    statement = '''
    Rscript %(R_SRC_PATH)s/bonferoni.R %(patterns)s --output %(outfile)s
    '''

    P.run(statement)

# }}}
# filter {{{
@transform(
    test_assoc,
    regex("^results/(.*)_assocs.tsv.gz$"),
    add_inputs(bonferoni),
    [r"results/\1_assocs_filtered.tsv", r"results/\1_stats.tsv"],
    r"\1"
    )
def filter(infiles, outfiles, pheno):

    assoc_gzip = infiles[0][0]
    bonf = infiles[1]
    filtered = outfiles[0]
    stats = outfiles[1]

    statement = '''
    echo "stat\\tvalue" > %(stats)s &&
    gzip -d -c %(assoc_gzip)s > %(pheno)s_temp.tsv &&
    wc -l %(pheno)s_temp.tsv | cut -f1 -d' ' | xargs -I @ echo -e 'kmers_tested\\t@' >> %(stats)s &&
    head -1 %(pheno)s_temp.tsv > %(filtered)s &&
    awk '$1=="bonf_thresh"{print $2}' %(bonf)s | xargs -I @ sh -c 'awk '\\''$4<@{print $0}'\\'' %(pheno)s_temp.tsv' >> %(filtered)s &&
    wc -l %(filtered)s | awk '{sum = $1 - 1; print sum}' | xargs -I @ echo -e 'significant_kmers\\t@' >> %(stats)s &&
    echo -e 'pheno\\t%(pheno)s' >> %(stats)s &&
    rm %(pheno)s_temp.tsv
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

    print(infile)

    ''' 
    Convert gff files into bed files. Also remove annotations that don't
    correspind to named genes (eg. hypothetical proteins) 
    '''

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

    

    statement = '''
    bwa index %(infile)s
    '''

    P.run(statement)

# # }}}

# map_kmers{{{
@transform(
    bonferoni,
    regex(r"^(.*)_stats.tsv$"),
    add_inputs(make_ref_list, ref2bed, annotation2bed, bwa_index),
    [r"\1_maps.tsv", r"\1_gene_info.tsv"],
    r"\1"
)
def map_kmers(infiles, outfiles, prefix):

    ''' Map kmers to annotation files '''

    to_cluster = False
    
    print(outfiles[0])
    print("\n#\n")
    print(outfiles[1])

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
    regex("^(.*)_maps.tsv$"),
    [r"\1_genes.tsv", r"\1_genes_near.tsv"],
    r"\1"
)
def summarise_genes(infiles, outfiles, prefix):

    ''' Count the number of kmer hits per gene, and calculate other
    statistics '''

    maps = infiles[0]
    gene_info = infiles[1]
    genes = outfiles[0]
    genes_near = outfiles[1]

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))

    statement = '''
    Rscript %(R_SRC_PATH)s/summarise_genes.R %(maps)s %(gene_info)s --prefix %(prefix)s &&
    awk '$1 != "significant_genes" && $1 != "significant_genes_near" {print $0}' %(prefix)s_stats.tsv > %(prefix)s_temp &&
    mv %(prefix)s_temp %(prefix)s_stats.tsv &&
    wc -l %(genes)s | awk '{sum = $1 - 1; print sum}' | xargs -I @ echo -e 'significant_genes\\t@' >> %(prefix)s_stats.tsv &&
    wc -l %(genes_near)s | awk '{sum = $1 - 1; print sum}' | xargs -I @ echo -e 'significant_genes_near\\t@' >> %(prefix)s_stats.tsv
    '''

    P.run(statement)

# }}}
# pathway_analysis {{{
@transform(
    summarise_genes,
    regex("^maps/(.*)_genes.tsv$"),
    r"maps/\1_pathways.tsv"
)
def pathwayAnalysis(infiles, outfile):
    pass

# }}}

# plot_genes {{{
@transform(
    summarise_genes,
    regex("^results/(.*)_genes.tsv$"),
    r"results/plots/\1.html",
    r"\1"
)
def plot_genes(infiles, outfile, pheno):

    template_path = os.path.dirname(os.path.realpath(__file__)) + "/plots/template.html"
    template = Template(filename=template_path)

    json_genes = ""
    with open(infiles[0]) as dat:
        titles = dat.readline().split(sep="\t")
        for l in dat:
            d = {}
            for t, f in zip(titles, l.split(sep="\t")):
                d[t.rstrip()] = f.rstrip()
            json_genes += json.dumps(d, indent=4) + ","
    page = template.render(pheno=pheno, json_genes=json_genes)
    html_file = open(outfile, "w")
    html_file.write(page)
    html_file.close()

    if not os.path.exists("results/plots/src"):
        shutil.copytree(os.path.dirname(os.path.realpath(__file__)) + "/plots/src", "results/plots/src")

# }}}
# summarise {{{
@follows(
    plot_genes
)
@merge(
    bonferoni,
    "results/summary.tsv"
)
def summarise(infiles, outfile):

    stats = [j for j in [i for sublist in infiles for i in sublist] if re.match(r"^.*stats.*$", j)]

    R_SRC_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "R"))
    
    statement = "Rscript %(R_SRC_PATH)s/summarise.R"

    for i in stats:
        statement += " " + i

    statement += " --output %(outfile)s"

    P.run(statement)

# }}}
# make_index {{{
@transform(
    summarise,
    regex(".*"),
    ["results/plots/index.html", "results/plots/src"]
    )
def make_index(infile, outfile):

    template_path = os.path.dirname(os.path.realpath(__file__)) + "/plots/index.html"
    template = Template(filename=template_path)

    table = "<table><thead><tr>"
    with open(infile, "r") as stats_file:
        header = stats_file.readline().rstrip().split("\t")
        table += "<th></th>"
        for h in header:
            table += "<th>" + h + "</th>"
        table += "</tr></thead><tbody>"
        for line in stats_file:
            stats = line.rstrip().split("\t")
            pheno = stats[0]
            table += "<tr><td><img src='" + pheno + "_density.png' style='width:20px;height:20px'></td><td><a href='" + pheno + ".html'>" + pheno + "</td>"
            for s in stats[1:]:
                table += "<td>" + s + "</td>"
            table += "</tr>"
        table += "</tbody></table>"

    page = template.render(table=table)
    html_file = open(outfile[0], "w")
    html_file.write(page)
    html_file.close()

    if os.path.exists(outfile[1]):
        shutil.rmtree(outfile[1])
    shutil.copytree(os.path.dirname(os.path.realpath(__file__)) + "/plots/src", outfile[1]);


# full {{{

@follows(
    summarise
)
def full():
    pass

# }}}

def main():
    P.main(sys.argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


