'''
Annotate Kmers
'''
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
    header = kmers_file.readline().rstrip()
    header = header + "\tgene_in\tgene_up\tgene_down\n"
    output_file.write(header)
    id_kmer = 0
    for line in kmers_file:
        id_kmer += 1
        kmers[id_kmer] = line.rstrip()

    kmers_file.close()
    print(str(len(kmers)) + " kmers loaded")

    # }}}

    for line in refs_file:

        query_fa_path = prefix + "_query.fa"
        query_bed_path = prefix + "_query.bed"
        query_sorted_bed_path = prefix + "_query_sorted.bed"

        if len(kmers) == 0:
            print("All kmers mapped!")
            break

        fa_path, bed_path = line.rstrip().split("\t")

        print("mapping to " + fa_path, end=": ")

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
            else:
                break

        # try to map kmers to reference index (bwa mem) {{{

        command = "bwa mem -v 1 -k 8 '" + fa_path + "' '" + query_fa_path + "'"
        results = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True) # stderr=subprocess.DEVNULL

        kmers_mapped = {}

        query_bed = open(query_bed_path, "w")

        for line in results.stdout:
            
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

            command = "bedtools intersect -a " + query_bed_path + " -b " + bed_path + " -wb"
            results_in = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
            genes_in = getGenes(results_in)

            query_sorted_bed = open(query_sorted_bed_path, "w")
            command = "bedtools sort -i " + query_bed_path + " > " + query_sorted_bed_path
            subprocess.run(command, shell=True, check=True)
            query_sorted_bed.close()
            
            command = "bedtools closest -a " + query_sorted_bed_path + " -b " + bed_path + " -D 'ref' -io -iu -nonamecheck"
            results_up = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
            genes_up = getGenes(results_up)

            command = "bedtools closest -a " + query_sorted_bed_path + " -b " + bed_path + " -D 'ref' -io -id -nonamecheck"
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

    if os.path.exists(query_fa_path):
        os.remove(query_fa_path)
    if os.path.exists(query_bed_path):
        os.remove(query_bed_path)
    if os.path.exists(query_sorted_bed_path):
        os.remove(query_sorted_bed_path)

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
