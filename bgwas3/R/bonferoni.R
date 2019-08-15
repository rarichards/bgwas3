# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(argparse)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(readr)));

description = """
Count unique kmer patterns and calculate p value threshold using bonferoni
Uses the patterns file generated from seer/pyseer
Requires unix 'sort'
"""
parser <- ArgumentParser(description = description);

parser$add_argument('file_patterns', help='patterns file');
parser$add_argument('out_file', help='out file (stats)');
parser$add_argument('--alpha', default=0.05, help='family-wise error rate');

args <- parser$parse_args()

statement <- paste0("LC_ALL=C sort -u -S 2014M ", args$file_patterns, " | wc -l");
unique_patterns <- strtoi(system(statement, intern=TRUE));

bonf_thresh <- args$alpha/unique_patterns;

data_frame(
		stat=c("unique_patterns", "bonf_thresh"),
		value=c(unique_patterns, bonf_thresh)
		) %>% 
	write_tsv(args$out_file);
