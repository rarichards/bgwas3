library(dplyr);
library(readr);
library(stringr);
library(tidyr);
library(ggtree);
library(tidytree);

"_static/tree.newick" %>% read.newick() -> dat_tree;

"phenos.tsv" %>% 
	read_tsv() -> dat;

dat %>% 
	select(id, clade_patient) %>%
	rename(
				 "Patient" = clade_patient
				 ) -> dat_temp;

rownames(dat_temp) <- dat %>% pull(id);

ggtree(dat_tree) -> p1;
p2 <- p1 %<+% dat_temp;

p2 + 
	geom_tiplab(size=2, align=TRUE, linesize=.5, offset=0.2) +
	geom_tippoint(aes(color=Patient)) -> p3;

normalize <- function(x) {
	return ((x - min(x)) / (max(x) - min(x)))
}

dat %>%
	select(id) %>%
	bind_cols(
						dat %>%
							select(starts_with("pheno"), -ends_with("int"), -ends_with("log")) %>%
							mutate_all(funs(normalize(.)))
	) %>%
	rename(
					"Tobromycin" = `pheno_tm`,
					"Impenem" = `pheno_ip`,
					"Aztreonam" = `pheno_at`,
					"Ciprofloxacin" = `pheno_ci`,
					"Colistin" = `pheno_co`,
					"Swim" = `pheno_swim`,
					"Swarm" = `pheno_swarm`,
					"Twitch" = `pheno_twitch`,
					"Hydrogen Cyanide" = `pheno_hcn`,
					"Cyanide" = `pheno_cyanide`,
					"2-Furoate" = `pheno_2furoate`,
					"3-Hydroxyisovalerate" = `pheno_3HIV`,
					"3-Methylthiopropionic acid" = `pheno_3MPA`,
					"Anthranilate" = `pheno_anthranilate`,
					"Betaine" = `pheno_betaine`,
					"Cystine" = `pheno_cystine`,
					"Formate" = `pheno_formate`,
					"Fumarate" = `pheno_fumarate`,
					"Histidine" = `pheno_histidine`,
					"Isoleucine" = `pheno_isoleucine`,
					"Leucine" = `pheno_leucine`,
					"Methanol" = `pheno_methanol`,
					"Methionine" = `pheno_methionine`,
					"Tryptophan" = `pheno_tryptophan`,
					"Uracil" = `pheno_uracil`,
					"Valine" = `pheno_valine`
	) -> dat;

rownames(dat) <- dat %>% pull(id);
dat %>% select(-id) -> dat;

gheatmap(
					p3,
					dat,
					offset = 0.3,
					hjust=0,
					width = 2,
					low = "white",
					high = "red",
					colnames = TRUE,
					colnames_position = "top",
					colnames_angle = 90,
					font.size = 3
) -> plot_tree;

png(filename="_static/tree.png", width=990, height=750);
plot_tree;
dev.off();
