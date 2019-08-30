library(dplyr);
library(readr);
library(reshape2);
library(stringr);
library(tidyr);
library(ggplot2);

"phenos.tsv" %>% 
	read_tsv() %>%
	select(starts_with("pheno"), -ends_with("int"), -ends_with("log")) %>%
	rename_all(funs(str_remove(., "pheno_"))) %>% 
	rename(
					"Tobromycin" = `tm`,
					"Impenem" = `ip`,
					"Aztreonam" = `at`,
					"Ciprofloxacin" = `ci`,
					"Colistin" = `co`,
					"Swim" = `swim`,
					"Swarm" = `swarm`,
					"Twitch" = `twitch`,
					"Hydrogen Cyanide" = `hcn`,
					"Cyanide" = `cyanide`,
					"2-Furoate" = `2furoate`,
					"3-Hydroxyisovalerate" = `3HIV`,
					"3-Methylthiopropionic acid" = `3MPA`,
					"Anthranilate" = `anthranilate`,
					"Betaine" = `betaine`,
					"Cystine" = `cystine`,
					"Formate" = `formate`,
					"Fumarate" = `fumarate`,
					"Histidine" = `histidine`,
					"Isoleucine" = `isoleucine`,
					"Leucine" = `leucine`,
					"Methanol" = `methanol`,
					"Methionine" = `methionine`,
					"Tryptophan" = `tryptophan`,
					"Uracil" = `uracil`,
					"Valine" = `valine`
) -> dat;

get_lower_tri<-function(cormat){
	cormat[upper.tri(cormat)] <- NA
	return(cormat)
}

dat %>%
	cor() %>%
	round(2) %>%
	get_lower_tri() %>%
	melt() %>%
	as_tibble() %>%
	rename(
				 "p1" = Var1,
				 "p2" = Var2,
				 "corr" = value
				 ) %>%
	ggplot(aes(x=p1, y=p2, fill=corr)) +
		geom_tile() +
		scale_fill_gradient2(
												 low = "blue", 
												 high = "red", 
												 mid = "white", 
												 midpoint = 0, 
												 limit = c(-1,1), 
												 space = "Lab", 
												 name="Pearson\nCorrelation"
		) +
		theme_minimal() + 
		theme(
					axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
					axis.title.x = element_blank(),
					axis.title.y = element_blank()
		) +
		coord_fixed() -> plot_cormat;
