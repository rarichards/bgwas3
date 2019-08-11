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
parser.add_argument("refs", help="text file listing annotation (paths of fa and bed)")
parser.add_argument("--prefix", help="prefix for output files")

args = parser.parse_args()
#args = parser.parse_args(["associations/log_at_assoc_filtered.txt", "ref.txt", "test.txt"])

# }}}

gene_info = {}

refs_file = open(args.refs, "r") 
output_file = open(args.prefix + "_maps.txt", "w")

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
        print("All kmers mapped!")
        break

    fa_path, bed_path = line.rstrip().split("\t")

    print("mapping to " + fa_path, end=": ")

    # command = "gff2bed < " + bed_path + " > reference.bed"
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
    results = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True, stderr=subprocess.DEVNULL)
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

                                if cigar_secondary == cigar:

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

    if kmers_mapped != False:

        def getGenes(results):
            genes_list = {}
            for line in results.stdout:
                fields = line.rstrip().split("\t")
                kmer_id, hit_id = fields[3].split("_")

                info = fields[9]
                gene = re.search("^.*gene=([^ ;]*);.*$", info)
                name = re.search("^.*name=([^ ;]*);.*$", info)
                description = re.search("^.*name=(\S* [^;]*);.*$", info)
                ID = re.search("^.*ID=([^ ;]*);.*$", info)

                if gene != None:
                    gene_name = gene.group(1)
                elif name != None:
                    gene_name = name.group(1)
                elif description != None:
                    gene_name = description.group(1)
                elif ID != None:
                    gene_name = ID.group(1)
                else:
                    break

                genes_list[kmer_id] = gene_name
                gene_info[gene_name] = info
            return genes_list

        command = "bedtools intersect -a query.bed -b " + bed_path + " -wb"
        results_in = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
        genes_in = getGenes(results_in)

        command = "bedtools sort -i query.bed > query_sorted.bed"
        subprocess.run(command, shell=True, check=True)
        
        command = 'bedtools closest -a query_sorted.bed -b ' + bed_path + ' -D "ref" -iu -nonamecheck'
        results_up = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
        genes_up = getGenes(results_up)

        command = 'bedtools closest -a query_sorted.bed -b ' + bed_path + ' -D "ref" -id -nonamecheck'
        results_down = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
        genes_down = getGenes(results_down)

        for id_kmer in kmers_mapped.keys():
            line = kmers_mapped[id_kmer] + "\t"
            if str(id_kmer) in genes_in:
                line += genes_in[str(id_kmer)]
            line += "\t"
            if str(id_kmer) in genes_up:
                line += genes_up[str(id_kmer)]
            line += "\t"
            if str(id_kmer) in genes_down:
                line += genes_down[str(id_kmer)]
            line += "\n"
            output_file.write(line)

    print(str(len(kmers_mapped)) + " new kmers mapped. " + str(len(kmers)) + " kmers left")

if os.path.exists("query.fa"):
    os.remove("query.fa")
if os.path.exists("query.bed"):
    os.remove("query.bed")
if os.path.exists("query_sorted.bed"):
    os.remove("query_sorted.bed")

output_file.close()
refs_file.close()   

with open(args.prefix + "_gene_info.txt", "w") as file_gene_info:
    for gene_name in gene_info.keys():
        line = gene_name + "\t" + gene_info[gene_name] + "\n"
        file_gene_info.write(line)

print("Finished mapping")
