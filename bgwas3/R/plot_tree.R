# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(readr)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(tidyr)));
suppressWarnings(suppressMessages(library(stringr)));
suppressWarnings(suppressMessages(library(argparse)));
suppressWarnings(suppressMessages(library(ggplot2)));
suppressWarnings(suppressMessages(library(ggtree)));
suppressWarnings(suppressMessages(library(tidytree)));
suppressWarnings(suppressMessages(library(ggnewscale)));
# suppressWarnings(suppressMessages(library(ggstance)));

get_args <- function(){

	suppressWarnings(suppressMessages(library(argparse)));

	description <- "make a range of phylogenetic tree and heatmap visualisations" 
	parser <- ArgumentParser(description = description);
	parser$add_argument("tree", help="newick file path");
	parser$add_argument("pheno", help="tsv file path");
	parser$add_argument("--output", help="output file path", default="tree.png");

 # return <- parser$parse_args()
 return(parser$parse_args(c("pangenome/accessory_binary_genes.fa.newick", "phenos.tsv")))
}

normalize <- function(x) {
	return ((x - min(x)) / (max(x) - min(x)))
}

main <- function(dat_pheno, dat_tree, output){

	dat_tree %>% 
		ggtree() +
		geom_tiplab(size=2, align=TRUE, linesize=.5) +
		geom_tippoint(aes=(color="patient")) -> p0;

	dat_pheno %>%
		select(id) %>%
		bind_cols(
							dat_pheno %>% 
								select(starts_with("clade")) %>%
								rename_all(funs(str_remove(., "clade_")))
		) -> dat_clades; 
	
	rownames(dat_clades) <- dat_clades %>% pull(id);

	dat_tree %>%
		as_tibble() %>%
		full_join(
			dat_clades %>% 
				rename("label" = id)
			)	%>% as.treedata() -> dat_tree2;

	dat_clades %>% select(-id) -> dat_clades;

	clade_name = colnames(dat_clades)[1]

	gheatmap(
						p0,
						dat_clades,
						offset = 0,
						width = 1,
						colnames = TRUE,
						colnames_position = "top",
						colnames_angle = 90,
						colnames_level = NULL,
						colnames_offset_x = 0,
						colnames_offset_y = 0,
						font.size = 4,
						hjust = 0.5
	) +
		scale_colour_brewer(palette = "Set1") -> p1;

palette = "Set1"
    scale_fill_viridis_d(option="patient", name="discrete\nvalue") -> p1;

	dat_pheno %>%
		select(id) %>%
		bind_cols(
							dat_pheno %>% 
								# select(starts_with("pheno"), -ends_with("log"), -ends_with("bin"), -ends_with("int")) %>%
								select(starts_with("pheno")) %>%
								mutate_all(funs(normalize(.))) %>%
								rename_all(funs(str_remove(., "pheno_")))
		) -> dat_phenos;

	rownames(dat_phenos) <- dat_phenos %>% pull(id);
	dat_phenos %>% select(-id) -> dat_phenos;

	p2 <- p1 + new_scale_fill();
	gheatmap(
						p2,
						dat_phenos,
						offset = 0,
						width = 1,
						low = "white",
						high = "red",
						colnames = TRUE,
						colnames_position = "top",
						colnames_angle = 90,
						colnames_level = NULL,
						colnames_offset_x = 0,
						colnames_offset_y = 0,
						font.size = 4,
						hjust = 0.5
	);
		# scale_fill_continuous();

	ggsave(plot=p, file=paste0(output));

}

args <- get_args();

args$pheno %>% read_tsv() -> dat_pheno;
args$tree %>% read.newick() -> dat_tree;

main(dat_pheno, dat_tree, args$output);
