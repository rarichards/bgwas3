# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(argparse, quietly=TRUE)));
suppressWarnings(suppressMessages(library(dplyr, quietly=TRUE)));
suppressWarnings(suppressMessages(library(readr, quietly=TRUE)));

description <- "Count the number of kmers mapped to each gene and calculate stats";
parser <- ArgumentParser(description = description); 
parser$add_argument("path_maps", help='maps file');
parser$add_argument('path_gene_info', help='gene info file');
parser$add_argument('path_output', help='output file');

args <- parser$parse_args()

args$path_maps %>% 
	read_tsv %>%
	group_by(gene_in) %>%
	mutate(
		"nlog10p" = -log10(`lrt-pvalue`),
		) %>%
	summarise(
		"hits" = n(),
		"min_beta" = min(beta),
		"max_beta" = max(beta),
		"mean_beta" = mean(beta),
		"var_beta" = var(beta),
		"min_nlog10p" = min(nlog10p),
		"max_nlog10p" = max(nlog10p),
		"mean_nlog10p" = mean(nlog10p),
		"var_nlog10p" = var(nlog10p),
		"min_af" = min(af),
		"max_af" = max(af),
		"mean_af" = mean(af),
		"var_af" = var(af)
		) %>%
	inner_join(
		args$path_gene_info %>% read_tsv(col_names=c("gene_in", "info"))
		) %>%
	write_tsv(args$path_output);
