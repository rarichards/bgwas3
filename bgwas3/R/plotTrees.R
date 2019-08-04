library(readr);
library(dplyr);
library(tidyr);
library(stringr);
library(argparse);
library(ggplot2);
library(ggtree);
library(tidytree);
library(ggstance);

parser <- ArgumentParser(description="Make a range of tree visualisations with a tree file and a phenotype file");
parser$add_argument("treefile", help="Tree file in newick format");
parser$add_argument("phenofile", help="Phenotype file in tsv format. Must include an 'id' column. Phenotypes must start with 'pheno_'. Clades must start with 'clade_'. Times must start with 'time_'");
parser$add_argument("outdir", help="Directory to output files", default=".");

args <- parser$parse_args();
#args <- parser$parse_args(c("tree.newick", "phenos.tsv", "out/static/trees"));

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

dat_tree %>% ggtree() -> p1;

if(length(clades) != 0){
	dat_tree %>% ggtree(aes(color=clades[0])) -> p1;
}

for(pheno in phenos){
	facet_plot(p1, panel=pheno, data=dat_phenos, geom=geom_barh, mapping=aes(x=dat_phenos[[pheno]]), stat="identity") -> p2;
	ggsave(file=paste0(args$outdir, "/", pheno, "tree.png"));
}
