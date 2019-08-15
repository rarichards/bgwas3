# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(argparse)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(readr)));
suppressWarnings(suppressMessages(library(ggplot2)));

description <- "draw a qqlot and histogram from pvalues";

parser <- ArgumentParser(description = description);

parser$add_argument('pvalues', help="list of p values");
parser$add_argument('--prefix', help="", default="plot");
parser$add_argument('--confidence-interval', help="", default=0.95);

args <- parser$parse_args();
#args <- parser$parse_args(c("temp1", "--confidence-interval=0.95"));

args$pvalues %>% 
	read_tsv(col_names=c("p"), col_types=cols(p=col_double()), skip=1) %>%
	arrange(p) %>%
	mutate(
		"observed" = -log10(p),
		"expected" = -log10(ppoints(n()))
		# "clower"   = -log10(qbeta(p = (1 - args$ci) / 2, shape1 = 1:n(), shape2 = n():1)),
    # "cupper"   = -log10(qbeta(p = (1 + args$ci) / 2, shape1 = 1:n(), shape2 = n():1))
		#"expected_nlog10p" = -log10(runif(n(), 0, 1))
		) -> dat;

dat %>% 
	ggplot(aes(x=expected, y=observed)) +
		geom_point() +
		geom_abline(intercept=0, slope=1);

dat %>% ggplot(aes(sample=p)) + stat_qq();

ggsave(paste0(args$prefix, "_qq.png"));

dat %>%
	ggplot(aes(x=p)) +
	geom_histogram();

ggsave(paste0(args$prefix, "_hist.png"));
