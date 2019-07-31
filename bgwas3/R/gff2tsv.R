library("dplyr");
library("tidyr");
library("readr");
library("stringr");

commandArgs(trailingOnly = TRUE) -> args;

args[1] -> file_in;
file_in %>% str_remove(".gff") %>% paste0(".tsv") -> file_out;

file_in %>% read_tsv(comment="#", col_names=c("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")) -> gff;

for(i in c("ID", "gene", "product", "db_xref", "inference", "protein_id")){
	regex <- paste0("(?<=", i, "=\")[^\"]*");
	gff %>%
		rowwise() %>%
		mutate(!!i := paste(unlist(str_extract_all(attributes, regex)), collapse=", ")) -> gff;
}

gff %>% 
	select(-c(attributes, source, score, strand, phase, feature, start, end, sequence)) %>%
	write_tsv(file_out)
