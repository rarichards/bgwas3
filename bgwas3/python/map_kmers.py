'''
Annotate Kmers '''
import sys
import os
import re
import subprocess
import argparse
import tempfile


def parse_aguments(): # {{{

    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("kmers_path", help="kmers file path (filtered output from SEER)")
    parser.add_argument("refs_path", help="text file listing annotation (paths of fa and bed)")
    parser.add_argument("--prefix", help="prefix for output files", default="out")

    return parser.parse_args()

# }}}

def main(kmers_path, refs_path, prefix):

    refs_file = open(refs_path, "r") 
    output_file = open(prefix + "_maps.tsv", "w")

    gene_info = {}
    kmers = {}

    # load kmers {{{
    kmers_file = open(kmers_path, "r")
    header = kmers_file.readline().rstrip("\n")
    header = header + "\tgene_in\tgene_up\tgene_down\n"
    output_file.write(header)
    id_kmer = 0
    for line in kmers_file:
        id_kmer += 1
        kmers[id_kmer] = line.rstrip("\n")

    kmers_file.close()
    print(str(len(kmers)) + " kmers loaded")

    # }}}

    for line in refs_file:

        query_fa_path = prefix + "_query.fa"
        query_bed_path = prefix + "_query.bed"
        query_sorted_bed_path = prefix + "_query_sorted.bed"
        mem_results_path = prefix + "_mem_results"

        if len(kmers) == 0:
            print("All kmers mapped!")
            break

        fa_path, bed_path = line.rstrip().split("\t")

        print("mapping to " + fa_path + ":")

        # make a fasta file with from unmapped kmers
        query_fa_path = prefix + "_query.fa"
        query_fa = open(query_fa_path, "w")
        for id_kmer in kmers.keys():
            query_fa.write(">" + str(id_kmer) + "\n" + kmers[id_kmer].rstrip().split("\t")[0] + "\n")
        query_fa.close()


        # bwa index
        for suffix in ["amb", "ann", "bwt", "pac", "sa"]:
            if not os.path.exists(fa_path + "." + suffix):
                command = "bwa index " + fa_path
                subprocess.run(command, shell=True, check=True)


        # try to map kmers to reference index (bwa mem) {{{

        if not os.path.exists(mem_results_path):
            command = "bwa mem " + fa_path + " " + query_fa_path + " > " + mem_results_path
            print(command)
            subprocess.run(command, shell=True, check=True)
            # results = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True) # stderr=subprocess.DEVNULL

        kmers_mapped = {}

        query_bed = open(query_bed_path, "w")
        mem_results = open(mem_results_path, "r")

        # for line in results.stdout:
        for line in mem_results:
            
            fields = line.rstrip().split("\t")

            if fields[0][0] == "@":
                continue

            if int(fields[1]) != 4:

                id_kmer = int(fields[0])

                kmers_mapped[id_kmer] = kmers[id_kmer]

                # del kmers[id_kmer]

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
                                        
                    except ValueError:
                        pass

                # }}}

        query_bed.close()

        for key in kmers_mapped.keys():
            del kmers[key] 

        # }}}

        if kmers_mapped != {}:

            def getGenes(results_path, info_index):
                genes_list = {}
                with open(results_path, encoding="utf8", errors='ignore') as results:
                    for line in results:
                    # for line in results.stdout:
                        fields = line.rstrip().split("\t")
                        kmer_id, hit_id = fields[3].split("_")

                        # info = fields[9]
                        info = fields[info_index]
                        # gene = re.search("^.*gene=([^ ;]*);.*$", info)
                        gene = re.findall(r'gene="([^ "]*)";', info)
                        name = re.findall(r'name="([^ "]*)";', info)
                        # name = re.search("^.*name=([^ ;]*);.*$", info)
                        description = re.search("^.*name=(\S* [^;]*);.*$", info)
                        ID = re.search("^.*ID=([^ ;]*);.*$", info)

                        if len(gene) != 0:
                            # gene_name = gene.group(1).strip('"')
                            gene_name = min(gene, key=len)
                        elif len(name) != 0:
                            # gene_name = name.group(1).strip('"')
                            gene_name = min(name, key=len)
                        elif description != None:
                            gene_name = description.group(1).strip('"')
                        elif ID != None:
                            gene_name = ID.group(1).strip('"')
                        else:
                            break

                        if kmer_id in genes_list.keys():
                            genes_list[kmer_id] += ";" + gene_name
                        else:
                            genes_list[kmer_id] = gene_name
                        gene_info[gene_name] = info

                # os.remove(results_path)

                return genes_list

            path_genes_in = prefix + "_genes_in"
            command = "bedtools intersect -a " + query_bed_path + " -b " + bed_path + " -wb > " + path_genes_in
            # results_in = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
            # genes_in = getGenes(results_in, -1)
            print(command)
            subprocess.run(command, shell=True)
            genes_in = getGenes(path_genes_in, -1)

            query_sorted_bed = open(query_sorted_bed_path, "w")
            command = "bedtools sort -i " + query_bed_path + " > " + query_sorted_bed_path
            subprocess.run(command, shell=True, check=True)
            query_sorted_bed.close()
            
            path_genes_up = prefix + "_genes_up"
            command = "bedtools closest -a " + query_sorted_bed_path + " -b " + bed_path + " -D 'ref' -io -iu -nonamecheck > " + path_genes_up
            print(command)
            # results_up = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
            subprocess.run(command, shell=True)
            genes_up = getGenes(path_genes_up, -2)

            path_genes_down = prefix + "_genes_down"
            command = "bedtools closest -a " + query_sorted_bed_path + " -b " + bed_path + " -D 'ref' -io -id -nonamecheck > " + path_genes_down
            print(command)
            # results_down = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
            subprocess.run(command, shell=True)
            genes_down = getGenes(path_genes_down, -2)

            print("genes_in: " + str(len(genes_in)))
            print("genes_up: " + str(len(genes_up)))
            print("genes_down: " + str(len(genes_down)))

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

        os.remove(mem_results_path);

    # if os.path.exists(query_fa_path):
    #     os.remove(query_fa_path)
    # if os.path.exists(query_bed_path):
    #     os.remove(query_bed_path)
    # if os.path.exists(query_sorted_bed_path):
    #     os.remove(query_sorted_bed_path)

    output_file.close()
    refs_file.close()   

    with open(prefix + "_gene_info.tsv", "w") as file_gene_info:
        for gene_name in gene_info.keys():
            line = gene_name + "\t" + gene_info[gene_name] + "\n"
            file_gene_info.write(line)

    print("Finished mapping")

# }}}

if __name__ == "__main__":
    args = parse_aguments()
    main(args.kmers_path, args.refs_path, args.prefix)
