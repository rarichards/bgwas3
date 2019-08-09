library(argparse);
library(dplyr);
library(readr);
library(ggplot2);

parser <- ArgumentParser(description="Draw a QQ-plot fro out of pyseer");
parser$add_argument('file_assoc', help='Associations file (output of pyseer in tsv format)');

#args <- parser$parse_args()
args <- parser$parse_args(c("temp.tsv"));

args$file_assoc %>% read_tsv() -> dat;

dat %>%
	mutate(
		"observed" = -log10(`lrt-pvalue`),
		"expected" = -log10(runif(n(), 0, 1))
		) %>%
	select(c("observed", "expected")) -> dat2;

dat2 %>% ggplot(aes(x=expected, y=observed)) +
	geom_point() +
	geom_abline(intercept=0, slope=1);
					
