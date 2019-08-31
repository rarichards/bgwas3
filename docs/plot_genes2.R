library(dplyr);
library(readr);
library(stringr);
library(tidyr);
library(ggplot2);
library(magick);

dat <- tibble(id = c(
												"3HIV", 
												"3HIV_log", 
												"anthranilate", 
												"anthranilate_log", "at_log", 
												"ci", 
												"co_log", 
												"isoleucine", 
												"isoleucine_log", 
												"leucine", 
												"leucine_log", 
												"methanol", 
												"methanol_log", 
												"swarm", 
												"swim_log", 
												"tm", 
												"tryptophan", 
												"tryptophan_log", 
												"uracil", 
												"uracil_log"
												)
) %>% mutate(name = "");

dat$name[dat$id == "3HIV"] <- "3-Hydroxyisovalerate";
dat$name[dat$id == "3HIV_log"] <- "log(3-Hydroxyisovalerate)";
dat$name[dat$id == "anthranilate"] <- "Anthranilate";
dat$name[dat$id == "anthranilate_log"] <- "log(Anthranilate)";
dat$name[dat$id == "at_log"] <- "log(Aztreonam)";
dat$name[dat$id == "ci"] <- "Ciprofloxacin";
dat$name[dat$id == "co_log"] <- "log(Colistin)";
dat$name[dat$id == "isoleucine"] <- "Isoleucine";
dat$name[dat$id == "isoleucine_log"] <- "log(Isoleucine)";
dat$name[dat$id == "leucine"] <- "Leucine";
dat$name[dat$id == "leucine_log"] <- "log(Leucine)";
dat$name[dat$id == "methanol"] <- "Methanol";
dat$name[dat$id == "methanol_log"] <- "log(Methanol)";
dat$name[dat$id == "swarm"] <- "Swarm";
dat$name[dat$id == "swim_log"] <- "log(Swim)";
dat$name[dat$id == "tm"] <- "Tobromycin";
dat$name[dat$id == "tryptophan"] <- "Tryptophan";
dat$name[dat$id == "tryptophan_log"] <- "log(Tryptophan)";
dat$name[dat$id == "uracil"] <- "Uracil";
dat$name[dat$id == "uracil_log"] <- "log(Uracil)";

dat_all <- tibble(
									Phenotype = as.character(),
									Gene = as.character(),
									Hits = as.numeric(),
									`Max -log10(p)` = as.numeric(),
									`Mean beta` = as.numeric(),
									`Mean MAF` = as.numeric()
									);

for(pheno in dat %>% pull(id)){
	paste0("_dat/", pheno, "_genes.tsv") %>% 
		read_tsv() %>%
		mutate(Phenotype = dat$name[dat$id == pheno])%>%
		rename(
					 "Gene" = gene,
					 "Mean beta" = mean_beta,
					 "Max -log10(p)" = max_nlog10p,
					 "Hits" = hits,
					 "Mean MAF" = mean_maf
		) %>%
		select(c("Phenotype", "Gene", "Hits", "Max -log10(p)", "Mean beta", "Mean MAF")) %>%
		ggplot(aes(x=`Mean beta`, y=`Max -log10(p)`, colour=`Mean MAF`, size=Hits, label=Gene)) +
		geom_point(alpha=0.5) +
		scale_size("Number of k-mers", range=c(1,10)) +
		scale_colour_gradient('Mean MAF') +
		theme_bw(base_size=14) +
		xlab("Mean effect size") +
		ylab("Maximum -log10(p-value)") +
		theme_classic() +
		theme(legend.position = "none") -> plot_gene;
	paste0("_static/", pheno, "_genes.png") %>%
		ggsave(plot, width=2, height=1.8, units="in", dpi=150, scale=2);
	plot_both <- image_append(image_scale(c(plot_gene, plot_qq), "500"));
	image_write(plot_both, paste0("_static/", pheno, "_both.png"));
}
