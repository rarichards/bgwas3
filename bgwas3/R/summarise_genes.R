# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(argparse, quietly=TRUE)));
suppressWarnings(suppressMessages(library(dplyr, quietly=TRUE)));
suppressWarnings(suppressMessages(library(readr, quietly=TRUE)));

description <- "Count the number of kmers mapped to each gene and calculate stats";
parser <- ArgumentParser(description = description); 
parser$add_argument("path_maps", help='maps file');
parser$add_argument('path_gene_info', help='gene info file');
parser$add_argument('path_output', help='output file');

args <- parser$parse_args();
#args <- parser$parse_args(c("results/swarm_mean_maps.txt", "results/swarm_mean_gene_info.txt", "test.tsv"));

tibble(`no genes mapped` = numeric()) -> dat_empty;

if(file.info(args$path_gene_info)$size == 0){
	dat_empty %>% write_tsv(args$path_output);
}else{

	args$path_maps %>% 
		read_tsv(col_types=cols(af="d", `filter-pvalue`="d", `lrt-pvalue`="d", beta="d", `beta-std-err`="d")) %>%
		group_by(gene_in) %>%
		mutate(
			"nlog10p" = -log10(`lrt-pvalue`),
			"maf" = if_else(af > 0.5, 1 - af, af)
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
			"min_maf" = min(maf),
			"max_maf" = max(maf),
			"mean_maf" = mean(maf),
			"var_maf" = var(maf)
			) %>%
		inner_join(
			args$path_gene_info %>% read_tsv(col_names=c("gene_in", "info"))
			) %>%
		write_tsv(args$path_output);
}
