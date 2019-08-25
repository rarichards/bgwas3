# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(argparse, quietly=TRUE)));
suppressWarnings(suppressMessages(library(dplyr, quietly=TRUE))); suppressWarnings(suppressMessages(library(readr, quietly=TRUE)));
suppressWarnings(suppressMessages(library(tidyr , quietly=TRUE)));

description <- "Count the number of kmers mapped to each gene and calculate stats";
parser <- ArgumentParser(description = description); 
parser$add_argument("path_maps", help='maps file');
parser$add_argument('path_gene_info', help='gene info file');
parser$add_argument('--prefix', help="prefix", default="test");

args <- parser$parse_args();
# args <- parser$parse_args(c("results/penicillin_maps.tsv", "results/penicillin_gene_info.tsv"));

tibble(`no genes mapped` = numeric()) -> dat_empty;

args$path_gene_info %>% read_tsv(col_names=c("gene", "info")) -> dat_gene_info;

args$path_maps %>% 
	read_tsv() %>%
	mutate(
		"nlog10p" = -log10(`lrt-pvalue`),
		"maf" = if_else(af > 0.5, 1 - af, af)
		) %>%
	rename(
				 "in" = gene_in,
				 "up" = gene_up,
				 "down" = gene_down
				 ) %>% 
	gather(key = "position", value = "gene", c("in", "up", "down")) %>%
	separate_rows(gene, sep=";") %>%
	filter(!is.na(gene)) -> dat_all;

dat_all %>% 
	filter(position == "in") %>%
	group_by(gene) %>%
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
						"min_maf" = min(maf),
						"max_maf" = max(maf),
						"mean_maf" = mean(maf),
						"var_maf" = var(maf)
						) %>%
	left_join(dat_gene_info) %>%
	arrange(desc(hits)) %>%
	write_tsv(paste0(args$prefix, "_genes.tsv"));

dat_all %>% 
	group_by(gene) %>%
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
						"min_maf" = min(maf),
						"max_maf" = max(maf),
						"mean_maf" = mean(maf),
						"var_maf" = var(maf)
						) %>%
	left_join(dat_gene_info) %>%
	arrange(desc(hits)) %>%
	write_tsv(paste0(args$prefix, "_genes_near.tsv"));
