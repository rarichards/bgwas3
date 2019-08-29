# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(argparse)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(readr)));

get_args <- function(){
	parser <- ArgumentParser();
	parser$add_argument('patterns', help='patterns file path', default="patterns.txt");
	parser$add_argument('--output', help='out file path', default="bonf.txt");
	parser$add_argument('--alpha', default=0.05, help='family-wise error rate');
	return(parser$parse_args());
}

main <- function(path_patterns, path_output, alpha){
	statement <- paste0("LC_ALL=C sort -u -S 2014M ", path_patterns, " | wc -l");
	unique_patterns <- strtoi(system(statement, intern=TRUE));
	bonf_thresh <- 0;
	if(unique_patterns != 0){
		bonf_thresh <- alpha/unique_patterns;
	}
	tibble(
			stat=c("unique_patterns", "bonf_thresh"),
			value=c(unique_patterns, bonf_thresh)
	) %>% 
		write_tsv(path_output);
}

args <- get_args();
path_patterns <- args$patterns;
path_output <- args$output;
alpha <- args$alpha;

main(args$patterns, args$output, args$alpha);
