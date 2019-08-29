# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(readr)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(tidyr)));
suppressWarnings(suppressMessages(library(stringr)));
suppressWarnings(suppressMessages(library(argparse)));
suppressWarnings(suppressMessages(library(ggplot2)));
suppressWarnings(suppressMessages(library(ggtree)));
suppressWarnings(suppressMessages(library(tidytree)));
suppressWarnings(suppressMessages(library(ggstance)));

description <- "Make a range of tree visualisations with a tree file and a phenotype file" 
parser <- ArgumentParser(description = description);
parser$add_argument("treefile", help="Tree file in newick format");
parser$add_argument("phenofile", help="Phenotype file in tsv format. Must include an 'id' column. Phenotypes must start with 'pheno_'. Clades must start with 'clade_'. Times must start with 'time_'");
parser$add_argument("outdir", help="Directory to output files", default=".");

# args <- parser$parse_args();
args <- parser$parse_args(c("pangenome/accessory_binary_genes.fa.newick", "phenos.tsv", "results/plot"));

message("# loading data #");
args$phenofile %>% read_tsv() -> dat_phenos;
dat_phenos %>% colnames() %>% str_subset("^pheno_.*") -> phenos;
dat_phenos %>% colnames() %>% str_subset("^clade_*") -> clades;
dat_phenos %>% colnames() %>% str_subset("^time_*") -> times;
args$treefile %>% read.newick() -> dat_tree;

dat_tree %>%
	as_tibble() %>%
 	full_join(
		dat_phenos %>% 
			rename("label" = id)
		)	%>% as.treedata() -> dat_tree;

message("# plotting trees #");
dat_tree %>% ggtree() -> p1;

if(length(clades) != 0){
	dat_tree %>% ggtree(aes(color=clades[0])) -> p1;
}

for(pheno in phenos){
	facet_plot(p1, panel=pheno, data=dat_phenos, geom=geom_barh, mapping=aes(x=dat_phenos[[pheno]]), stat="identity") -> p2;
	ggsave(file=paste0(args$outdir, "/", pheno, "tree.png"));
}

"pangenome/accessory_binary_genes.fa.newick" %>% read.newick -> tree;
tree %>% ggtree() -> p1;

normalize <- function(x) {
return ((x - min(x)) / (max(x) - min(x)))
}

"phenos.tsv" %>% 
	read_tsv %>% 
	select(id) %>%
	bind_cols(
						"phenos.tsv" %>% 
							read_tsv() %>%
							select(starts_with("pheno"), -ends_with("log"), -ends_with("bin"), -ends_with("int")) %>%
							mutate_all(funs(normalize(.))) %>%
							rename_all(funs(str_remove(., "pheno_")))
	) -> dat; 

rownames(dat) <- dat %>% pull(id);

dat %>% select(-id) -> dat2;


gheatmap(p1, dat2, offset = 0, width = 1, low = "green",
high = "red", color = "white", colnames = TRUE,
colnames_position = "bottom", colnames_angle = 90,
colnames_level = NULL, colnames_offset_x = 0,
colnames_offset_y = 0, font.size = 4, hjust = 0.5)





	gather(key="key", value="value", -id) %>%
	ggplot(aes(x=key, y=id, fill=value)) + geom_tile() -> p2;


	ggplot(aes(y=id, fill=pheno_at)) + geom_tile();
