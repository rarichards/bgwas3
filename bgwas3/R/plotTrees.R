suppressPackageStartupMessages(library(readr));
suppressPackageStartupMessages(library(dplyr));
suppressPackageStartupMessages(library(tidyr));
suppressPackageStartupMessages(library(stringr));
suppressPackageStartupMessages(library(argparse));
suppressPackageStartupMessages(library(ggplot2));
suppressPackageStartupMessages(library(ggtree));
suppressPackageStartupMessages(library(tidytree));

parser <- ArgumentParser(description="Make a range of tree visualisations with a tree file and a phenotype file");
parser$add_argument("treefile", help="Tree file in newick format");
parser$add_argument("phenofile", help="Phenotype file in tsv format. Must include an 'id' column. Phenotypes must start with 'pheno_'. Clades must start with 'clade_'. Times must start with 'time_'");
parser$add_argument("outdir", help="Directory to output files", default=".");

#args <- parser$parse_args();
args <- parser$parse_args(c("tree.newick", "phenos.tsv", "plots"));

args$phenofile %>% read_tsv() -> dat_phenos;
dat_phenos %>% colnames() %>% str_subset("^pheno_.*") -> phenos;
dat_phenos %>% colnames() %>% str_subset("^clade_*") -> clades;
dat_phenos %>% colnames() %>% str_subset("^time_*") -> times;

args$treefile %>%
	read.newick() %>%
	as_tibble() %>%
 	full_join(
		dat_phenos %>% 
			rename("label" = id)
		) %>% as.treedata() -> dat_tree;

if(length(clades) == 0){
	for(pheno in phenos){
		dat_tree %>%
			ggtree() +
				geom_tiplab() +
				geom_facet(panel = pheno, data = dat_phenos, 
					geom = ggstance::geom_barh, 
          aes(x = pheno), 
          stat = "identity", width = .6
				) + 
				theme_tree2(legend.position=c(.05, .85));
			# ggsave(paste0(args$outdir, "/", clade, pheno, "tree.png"));
	}


for(clades in clades){
	for(pheno in phenos){
		dat_tree %>%
			ggtree(aes(color=pheno)) +
				geom_tiplab() +
				geom_facet(panel = pheno, data = dat_phenos, 
					geom = ggstance::geom_barh, 
          aes(x = pheno, color = clade, fill = clade), 
          stat = "identity", width = .6
				) + 
				theme_tree2(legend.position=c(.05, .85));
			# ggsave(paste0(args$outdir, "/", clade, pheno, "tree.png"));
	}
}

dat_tree %>%
	ggtree(aes(color=`pheno_penicillin`)) +
		geom_tippoint(aes(color=`pheno_penicillin`)) +
		geom_tiplab();
