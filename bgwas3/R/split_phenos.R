# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(readr)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(tidyr)));
suppressWarnings(suppressMessages(library(argparse)));
suppressWarnings(suppressMessages(library(stringr)));

parser <- ArgumentParser(description="Extract phenotype columns from a tsv file and write seperate tsv files for each");
parser$add_argument("infile", help="File in tsv format that includes an 'id' column and one or more 'pheno_*' columns");
parser$add_argument("outdir", help="Directory to write new tsv files", default=".");

args <- parser$parse_args();

args$infile %>% 
	read_tsv() -> dat;

for(col in dat %>% select(starts_with("pheno_")) %>% colnames()){
	pheno = str_sub(col, 7, -1);
	paste(pheno);
	dat %>% 
		select(c("id", col)) %>%
		rename(pheno = col) %>%
		write_tsv(paste0(args$outdir, "/", pheno, ".tsv"));
}
