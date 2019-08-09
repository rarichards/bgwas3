
# bwaIndex {{{
@transform(
    assembly,
    suffix(".fa"),
    [".amb", ".ann", ".bwt", ".pac", ".sa"]
    )
def bwaIndex(infile, outfiles):
    
    statement = '''
    bwa index %(infile)s
    '''

    P.run(statement)

# }}}
# makeRefList {{{
@follows(
    mkdir("refs")
    )
@merge(
    [annotateGenomes, [assembly], ["refs/*"]],
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

# seer2fa {{{
@transform(
    bonferoniFilter,
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
    [testAssociations, bonferoniFilter],
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

