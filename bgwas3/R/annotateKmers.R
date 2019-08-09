library(readr);
library(dplyr);
library(tidyr);
library(argparse);
library(ggplot2);

parser <- ArgumentParser(description="");
parser$add_argument("kmers", help="");
parser$add_argument("references", help="");
parser$add_argument("output", help="");

args <- parser$parse_args();

args$references %>% read_tsv() -> dat_references;

args$kmers %>% file(open="r") -> file_kmers;
"remaining_kmers.txt" %>% file(open="w") -> file_remaining_kmers_fa;

count <- 0
while (length(line <- readLines(file_kmers, n = 1, warn = FALSE)) > 0) {
	strsplit(line, "\t")[1] -> fa;
	count <- count + 1;
	pasete0(">", count, "\n", fa, "\n")
} 



