suppressWarnings(suppressMessages(library(argparse)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(readr)));
suppressWarnings(suppressMessages(library(tidyr)));

parser <- ArgumentParser();
parser$add_argument('tsv', nargs="+");
parser$add_argument('--output', help="", default="summary.tsv")

args <- parser$parse_args();

dat <- tibble(
	"pheno" = character(),
	"bonf_thresh" = numeric(),
	"significant_kmers" = numeric(),
	"significant_genes" = numeric()
)

for(tsv in args$tsv){
	dat %>%
		bind_rows(
			tsv %>%
				read_tsv() %>%
				spread(stat, value) %>%
				select(c("pheno", "bonf_thresh", "significant_kmers", "significant_genes")) %>%
				mutate(
					"bonf_thresh" = as.numeric(bonf_thresh),
					"significant_kmers" = as.numeric(significant_kmers),
					"significant_genes" = as.numeric(significant_genes)
					)
			) -> dat;
}

write_tsv(dat, args$output);
