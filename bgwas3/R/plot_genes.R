# 2019 Gregory Leeman g-r-eg@outlook.com

suppressWarnings(suppressMessages(library(readr)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(tidyr)));
suppressWarnings(suppressMessages(library(ggplot2)));

commandArgs(trailingOnly = TRUE) -> args;

hits = args[1]
png = args[2]

hits %>% 
	read_tsv() %>%
	ggplot(aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene)) +
		 geom_point(alpha=0.5) +
		 #geom_text_repel(aes(size=60), show.legend = FALSE, colour='black') +
		 geom_text(aes(size=60), show.legend = FALSE, colour='black') +
		 scale_size("Number of k-mers", range=c(1,10)) +
		 scale_colour_gradient('Average MAF') +
		 theme_bw(base_size=14) +
		 ggtitle("Penicillin resistance") +
		 xlab("Average effect size") +
		 ylab("Maximum -log10(p-value)");
	
ggsave(png)
