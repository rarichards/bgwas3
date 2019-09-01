library(dplyr);
library(readr);
library(stringr);
library(tidyr);
library(ggplot2);
library(magick);

dat <- tibble(id = c(
												"3HIV", 
												"anthranilate", 
												"at_log", 
												"ci", 
												"co_log", 
												"isoleucine", 
												"leucine", 
												"methanol", 
												"swarm", 
												"swim_log", 
												"tm", 
												"tryptophan", 
												"uracil"
												)
) %>% mutate(name = "");

dat$name[dat$id == "3HIV"] <- "3-Hydroxyisovalerate";
dat$name[dat$id == "anthranilate"] <- "Anthranilate";
dat$name[dat$id == "at_log"] <- "Aztreonam";
dat$name[dat$id == "ci"] <- "Ciprofloxacin";
dat$name[dat$id == "co_log"] <- "Colistin";
dat$name[dat$id == "isoleucine"] <- "Isoleucine";
dat$name[dat$id == "leucine"] <- "Leucine";
dat$name[dat$id == "methanol"] <- "Methanol";
dat$name[dat$id == "swarm"] <- "Swarm";
dat$name[dat$id == "swim_log"] <- "Swim";
dat$name[dat$id == "tm"] <- "Tobromycin";
dat$name[dat$id == "tryptophan"] <- "Tryptophan";
dat$name[dat$id == "uracil"] <- "Uracil";

dat_all <- tibble(
									Phenotype = as.character(),
									Gene = as.character(),
									Hits = as.numeric(),
									`Max -log10(p)` = as.numeric(),
									`Mean beta` = as.numeric(),
									`Mean MAF` = as.numeric()
									);

p <- list();
i = 0;
for(pheno in dat %>% pull(id)){
	i = i + 1;
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
		ggsave(plot_gene, width=2, height=1.8, units="in", dpi=150, scale=2);

	image_read(paste0("_static/", pheno, "_genes.png")) -> plot_gene;
	image_read(paste0("_static/qq/", pheno, "_qq.png")) -> plot_qq;

	image_append(image_scale(c(plot_gene, plot_qq), "600")) %>%
		image_border("white", "50x50") %>%
		image_annotate(paste0(dat$name[dat$id == pheno]), font = 'Times', size = 50, location = "+100+0") -> plot_both;
	# image_write(plot_both, paste0("_static/", pheno, "_both.png"));
	# plot_both <- image_read(paste0("_static/", pheno, "both.png"));
	p[[i]] <- plot_both;
}

p[[14]] <- image_read("_static/legend.png");

plot_1 <- image_append(c(
														p[[1]],
														p[[2]],
														p[[3]],
														p[[4]],
														p[[5]],
														p[[6]],
														p[[7]]
														), stack=TRUE);

plot_2 <- image_append(c(
														p[[8]],
														p[[9]],
														p[[10]],
														p[[11]],
														p[[12]],
														p[[13]],
														p[[14]]
														), stack=TRUE);

plot_all <- image_append(c(plot_1, plot_2));

# plot_1 <- image_append(c(
# 														p[[1]],
# 														p[[2]],
# 														p[[3]],
# 														p[[4]],
# 														p[[5]]
# 														), stack=TRUE);
# plot_2 <- image_append(c(
# 														p[[6]],
# 														p[[7]],
# 														p[[8]],
# 														p[[9]],
# 														p[[10]]
# 														), stack=TRUE);
# plot_3 <- image_append(c(
# 														p[[11]],
# 														p[[12]],
# 														p[[13]],
# 														p[[14]],
# 														p[[15]]
# 														), stack=TRUE);
# plot_4 <- image_append(c(
# 														p[[16]],
# 														p[[17]],
# 														p[[18]],
# 														p[[19]],
# 														p[[20]]
# 														), stack=TRUE);
# plot_all <- image_append(c(plot_1, plot_2, plot_3, plot_4));

image_write(plot_all , "all.png");
