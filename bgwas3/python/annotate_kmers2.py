import sys
import os
import re
import tempfile
import subprocess
import pybedtools
import argparse

parser = argparse.ArgumentParser(description='Annotate Kmers')

parser.add_argument("kmers_file", help="Kmers file, filtered output from SEER")
parser.add_argument("gff", help="gff reference file")
parser.add_argument("output_file", help="output file")
parser.add_argument("--type" defult = "ref")

args = parser.parse_args()

output_file = open(options.output, "rw")

kmers_file = open(args.kmers, "r")
header = kmers_file.readline()
kmers_sum = 0
kmers_fa_file = open("remaining_kmers.fa", "w")
for kmer in kmers_file:
    kmers_sum += 1
    kmers_fa_file.write(">" + str(kmer_sum) + "\n" + kmer.split("\t")[0] _ "\t")

def bwa_mem(ref_fa, query_fa): # {{{
    command = "bwa mem -v 1 -k 8 '" + reference + "' '" + fasta + "'"
    bwa_p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    for sam_line in bwa_p.stdout:
        sam_fields = sam_line.rstrip().split("\t")

        # discard header
        if sam_fields[0][0] == "@":
            continue

        positions = []
        if int(sam_fields[1]) & 4 == 4:
            mapped = False
        else:
            mapped = True

            # primary mapping
            if int(sam_fields[1]) & 16 == 16:
                strand = "-"
            else:
                strand = "+"
            if len(sam_fields) < 10:
                mapped = False
                positions = True
            else:
                positions.append((sam_fields[2], sam_fields[3], int(sam_fields[3]) + len(sam_fields[9]) - 1, strand))

                # secondary mappings (as good as primary - same CIGAR string)
                if len(sam_fields) > 15:
                    try:
                        secondary = sam_fields[15].split(":")
                        if secondary[0] == "XA" and secondary[1] == "Z":
                            for secondary_mapping in secondary[2].split(";"):
                                if secondary_mapping != '':
                                    (contig, pos, cigar, edit_distance) = secondary_mapping.split(",")
                                    if cigar == sam_fields[5]:
                                        strand = pos[0]
                                        positions.append((contig, pos[1:], int(pos[1:]) + len(sam_fields[9]) - 1, strand))
                    # Ignore secondary mappings which don't match the expected format
                    except ValueError:
                        pass
        yield(mapped, positions)

def bwa_fastmap(ref_fa, query_fa): # {{{
    else:
        mapped = False
        positions = []

        first_line = bwa_p.stdout.readline().rstrip().split("\t")
        if first_line == ['']:
            return
        (sq, idx, length) = first_line
        while True:
            fastmap_line = bwa_p.stdout.readline()
            fastmap_line = fastmap_line.rstrip()
            if fastmap_line == "//":
                next_line = bwa_p.stdout.readline().rstrip().split("\t")
                fastmap_hit = BWA(mapped, positions)
                if len(next_line) < 3:  # EOF reached
                    yield(fastmap_hit)
                    return
                else:
                    (sq, idx, length) = next_line
                    mapped = False
                    positions = []
                    yield(fastmap_hit)
            else:
                hits = []
                fastmap_fields = fastmap_line.split("\t")
                # in case a line is missing a few fields
                if len(fastmap_fields) < 5 or fastmap_fields[4] == '*':
                    continue
                #
                if fastmap_fields[1] == '0' and fastmap_fields[2] == length: #  full hits only
                    mapped = True
                    for hit in fastmap_fields[4:]:
                        try:
                            (contig, pos) = hit.split(":")
                        except:
                            print(fastmap_fields[4:])
                        strand = pos[0]
                        positions.append((contig, int(pos[1:]), int(pos[1:]) + int(length) - 1, strand))





