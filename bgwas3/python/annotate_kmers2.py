import sys
import os
import re
import tempfile
import subprocess
import pybedtools 
import argparse

# argparse {{{
parser = argparse.ArgumentParser(description='Annotate Kmers')

parser.add_argument("kmers", help="kmers file path (filtered output from SEER)")
parser.add_argument("refs", help="text file listing annotation (paths of fa and gff)")
parser.add_argument("output", help="output file")

args = parser.parse_args()
#args = parser.parse_args(["associations/log_at_assoc_filtered.txt", "ref.txt", "test.txt"])

# }}}

refs_file = open(args.refs, "r")
output_file = open(args.output, "w")

# load kmers {{{
kmers_file = open(args.kmers, "r")
header = kmers_file.readline().rstrip()
header = header + "\tgene_in\tgene_up\tgene_down\n"
output_file.write(header)
kmers = {}
id_kmer = 0
for line in kmers_file:
    id_kmer += 1
    kmers[id_kmer] = line.rstrip()

kmers_file.close()
print(str(len(kmers)) + " kmers loaded")

# }}}

for line in refs_file:

    if len(kmers) == 0:
        print("all kmers mapped")
        break

    fa_path, gff_path = line.rstrip().split("\t")

    print("mapping to " + fa_path)
    # make index files from reference fasta (bwa index) (required for bwa mem)
    # command = "bwa index " + fa_path
    # subprocess.run(command, shell=True, check=True)

    # make a fasta file with from unmapped kmers
    query_fa = open("query.fa", "w")
    for id_kmer in kmers.keys():
        query_fa.write(">" + str(id_kmer) + "\n" + kmers[id_kmer].rstrip().split("\t")[0] + "\n")
        # if len(kmers[id_kmer]["maps"]) == 0:
        #     kmers_file_fa.write(">" + str(key) + "\n" + kmers[key]["seq"] + "\n")
    query_fa.close()

    query_bed = open("query.bed", "w")

    # try to map kmers to reference index (bwa mem) {{{
    command = "bwa mem -v 1 -k 8 '" + fa_path + "' 'query.fa'"
    results = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    # command2 = "bwa mem -v 1 -k 8 '" + fa_path + "' 'query.fa' > test"
    # subprocess.run(command2, shell=True, check=True)

    kmers_mapped = {}

    for line in results.stdout:
        fields = line.rstrip().split("\t")

        if fields[0][0] == "@":
            continue

        if int(fields[1]) != 4:

            id_kmer = int(fields[0])
            kmers_mapped[id_kmer] = kmers[id_kmer]
            del kmers[id_kmer]

            # primary mapping {{{
            id_hit = 1
            if int(fields[1]) == 16:
                strand = "-"
            else:
                strand = "+"
            contig = fields[2]
            start = int(fields[3])
            length = len(fields[9])
            end = start + length -1
            cigar = fields[5]

            # kmers[id_kmer]["maps"].append({})
            # kmers[id_kmer]["maps"][-1]["contig"] = contig
            # kmers[id_kmer]["maps"][-1]["start"] = start
            # kmers[id_kmer]["maps"][-1]["length"] = length
            # kmers[id_kmer]["maps"][-1]["end"] = end

            query_bed.write('\t'.join([contig, str(start), str(end), str(id_kmer) + "_" + str(id_hit), '0', strand]) + "\n")

            # }}}

            # secondary mappings {{{
            if len(fields) > 15:

                try:
                    secondary = fields[15].split(":")

                    if secondary[0] == "XA" and secondary[1] == "Z":

                        mappings = secondary[2].split(";")

                        for mapping in mappings:

                            if mapping != '':

                                id_hit += 1
                                fields = mapping.split(",")
                                strand = fields[1][0]
                                contig = fields[0]
                                start = int(fields[1][1:])
                                end = start + length -1
                                cigar_secondary = fields[2]

                                if secondary_cigar == cigar:

                                    query_bed.write('\t'.join([contig, str(start), str(end), str(id_kmer) + "_" + str(id_hit), '0', strand]) + "\n")

                                    # kmers[id_kmer]["maps"].append({})
                                    # kmers[id_kmer]["maps"][-1]["contig"] = contig
                                    # kmers[id_kmer]["maps"][-1]["start"] = start
                                    # kmers[id_kmer]["maps"][-1]["length"] = length
                                    # kmers[id_kmer]["maps"][-1]["end"] = end

                except ValueError:
                    pass

            # }}}


    # }}}

    query_bed.close()
    command = "bedtools sort -i query.bed > query_sorted.bed"
    subprocess.run(command, shell=True, check=True)

    genes_in = {}
    command = "bedtools intersect -a query.bed -b " + gff_path + " -wb"
    print(command)
    results = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    for line in results.stdout:
        fields = line.rstrip().split("\t")
        kmer_id, hit_id = fields[3].split("_")
        features = dict(map(lambda s : s.split('='), fields[14].split(";")))
        genes_in[kmer_id] = features["gene"]

    # genes_up = {}
    # command = 'bedtools closest -a query_sorted.bed -b ' + gff_path + ' -D "ref" -iu -nonamecheck'
    # print(command)
    # results = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    # for line in results.stdout:
    #     fields = line.rstrip().split("\t")
    #     kmer_id, hit_id = fields[3].split("_")
    #     features = dict(map(lambda s : s.split('='), fields[14].split(";")))
    #     genes_up[kmer_id] = features["gene"]

    # genes_down = {}
    # command = 'bedtools closest -a query_sorted.bed -b ' + gff_path + ' -D "ref" -id -nonamecheck'
    # results = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    # for line in results.stdout:
    #     fields = line.rstrip().split("\t")
    #     kmer_id, hit_id = fields[3].split("_")
    #     features = dict(map(lambda s : s.split('='), fields[14].split(";")))
    #     genes_down[kmer_id] = features["gene"]

    print(genes_in)
    for id_kmer in kmers_mapped.keys():
        line = kmers_mapped[id_kmer] + "\t"
        if str(id_kmer) in genes_in:
            line += genes_in[str(id_kmer)]
        # line += "\t"
        # if id_kmer in genes_up:
        #     line += genes_up[id_kmer]
        # line += "\t"
        # if id_kmer in genes_down:
        #     line += genes_down[id_kmer]
        line += "\n"
        output_file.write(line)

    print(str(len(kmers)) + " kmers left")
output_file.close()
refs_file.close()
