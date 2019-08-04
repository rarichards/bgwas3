suppressPackageStartupMessages(library(readr));
suppressPackageStartupMessages(library(dplyr));
suppressPackageStartupMessages(library(tidyr));
suppressPackageStartupMessages(library(argparse));
suppressPackageStartupMessages(library(stringr));

parser <- ArgumentParser(description="Extract phenotype columns from a tsv file and write seperate tsv files for each");
parser$add_argument("infile", help="File in tsv format that includes an 'id' column and one or more 'pheno_*' columns");
parser$add_argument("outdir", help="Directory to write new tsv files", default=".");

args <- parser$parse_args();

args$infile %>% 
	read_tsv() -> dat;

for(col in dat %>% select(starts_with("pheno_")) %>% colnames()){
	pheno = str_sub(col, 7, -1);
	dat %>% 
		select(c("id", col)) %>%
		rename(pheno = col) %>%
		write_tsv(paste0(args$outdir, "/", pheno, ".tsv"));
}
