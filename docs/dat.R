library(dplyr);
library(readr);
library(stringr);
library(tidyr);
library(reshape2);
library(ggtree);
library(tidytree);
library(ggplot2);


# tables {{{
dat_non_met = tibble(
									 Name = c(
														"Tobromycin",
														"Impenem",
														"Aztreonam",
														"Ciprofloxacin",
														"Colistin",
														"Swim",
														"Swarm",
														"Twitch"
									 )
) %>% mutate("Description" = "");

dat_met = tibble(
									 Chemical = c(
														"Hydrogen Cyanide",
														"Cyanide",
														"2-Furoate",
														"3-Hydroxyisovalerate",
														"3-Methylthiopropionic acid",
														"Anthranilate",
														"Betaine",
														"Cystine",
														"Formate",
														"Fumarate",
														"Histidine",
														"Isoleucine",
														"Leucine",
														"Methanol",
														"Methionine",
														"Tryptophan",
														"Uracil",
														"Valine"
									 )
);

dat_non_met$Description[dat_non_met$Name == "Tobromycin"] <- "Resistance to inhalant antibiotic Tobromycin" 
dat_non_met$Description[dat_non_met$Name == "Impenem"] <- "Resistance to intravenous antibiotic Imipenem";
dat_non_met$Description[dat_non_met$Name == "Aztreonam"] <- "Resistance to intravenous/intramuscular antibiotic Aztreonam";
dat_non_met$Description[dat_non_met$Name == "Ciprofloxacin"] <- "Resistance to oral antibiotic Ciprofloxacin";
dat_non_met$Description[dat_non_met$Name == "Colistin"] <- "Resistance to 'last-resort' antibioti Colistin";
dat_non_met$Description[dat_non_met$Name == "Swim"] <- "Measure of cell surface bactera movement by flagella";
dat_non_met$Description[dat_non_met$Name == "Swarm"] <- "Mesaure of rapid surface movement by multiple bacteria with rotating flagella";
dat_non_met$Description[dat_non_met$Name == "Twitch"] <- "Measure of slow baceria movement powered by pili";
dat_non_met$Description[dat_non_met$Name == "Hydrogen Cyanide"] <- "P. aeruginosa is one of a limited number of organisms which can synthesise cyanide, though the biological reason is unclear";

# }}}
# cormat {{{
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

# }}}
# density plots {{{
"phenos.tsv" %>% 
	read_tsv() %>%
	select(starts_with("pheno")) -> dat;

multiplot <- function(plots, cols=1, byrow=TRUE) { # {{{
	library(grid)

	length(plots) -> n_plots;

	if (n_plots==1) {
    print(plots[[1]])
		return;
	}

	matrix(
				 seq(1, cols * ceiling(n_plots/cols)),
				 ncol = cols,
				 nrow = ceiling(n_plots/cols),
				 byrow = byrow
	) -> layout;

	grid.newpage()
	grid::pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

	for (i in 1:n_plots) {
		as.data.frame(which(layout == i, arr.ind = TRUE)) -> matchidx;
		print(plots[[i]], vp = viewport(
																		layout.pos.row = matchidx$row,
																		layout.pos.col = matchidx$col
																		)
		);
  }


} # }}}

phenos <- c(
"pheno_at",
"pheno_tm",
"pheno_ip",
"pheno_ci",
"pheno_co",
"pheno_swim",
"pheno_swarm",
"pheno_twitch",
"pheno_hcn",
"pheno_cyanide",
"pheno_2furoate",
"pheno_3HIV",
"pheno_3MPA",
"pheno_anthranilate",
"pheno_betaine",
"pheno_cystine",
"pheno_formate",
"pheno_fumarate",
"pheno_histidine",
"pheno_isoleucine",
"pheno_leucine",
"pheno_methanol",
"pheno_methionine",
"pheno_tryptophan",
"pheno_uracil",
"pheno_valine"
);

i <- 0;
p <- list();

for(pheno in phenos){ # {{{

	log <- paste0(pheno, "_log");
	int <- paste0(pheno, "_int");
	pheno %>% 
		str_remove("pheno_") %>% 
		str_to_title() -> name;

	if(name == "3hiv"){
		name <- "3HIV";
	}
	if(name == "3Mpa"){
		name <- "3MPA";
	}

	dat_temp <- tibble(
											pheno = dat %>% pull(pheno),
											log = dat %>% pull(log),
											int = dat %>% pull(int)
											);

	i <- i + 1;
	bins = nclass.FD(dat_temp[["pheno"]]);
	p[[i]] <- eval(
						 dat_temp %>% 
							 ggplot(aes(x=pheno)) +
								geom_density(fill="black") +
								geom_vline(aes(xintercept=mean(pheno))) +
								# geom_histogram(bins=bins, fill="black") +
								xlab(name) +
								# ylab("Count") +
								theme_classic() +
								theme(
									axis.title.y=element_blank()
								)
						);

	i <- i + 1;
	bins = nclass.FD(dat_temp[["log"]]);
	p[[i]] <- eval(
						 dat_temp %>% 
							 ggplot(aes(x=log)) +
								geom_density(fill="black") +
								geom_vline(aes(xintercept=mean(log))) +
								# geom_histogram(bins=bins, fill="black") +
								xlab(paste0("log10(", name, ")")) +
								theme_classic() +
								theme(
									axis.title.y=element_blank()
								)
						);

	i <- i + 1;
	bins = nclass.FD(dat_temp[["int"]]);
	p[[i]] <- eval(
						 dat_temp %>% 
							 ggplot(aes(x=int)) +
								geom_density(fill="black") +
								geom_vline(aes(xintercept=mean(int))) +
								# geom_histogram(bins=bins, fill="black") +
								xlab(paste0("INT(", name, ")")) +
								theme_classic() +
								theme(
									axis.title.y=element_blank()
								)
						);
}

png(filename="_static/dens.png", width=750, height=990);
p %>% multiplot(cols = 9) -> plot_dens;
dev.off();

# }}}

# for(pheno in phenos){ # {{{

# 	log <- paste0(pheno, "_log");
# 	int <- paste0(pheno, "_int");
# 	pheno %>% 
# 		str_remove("pheno_") %>% 
# 		str_to_title() -> name;

# 	if(name == "3hiv"){
# 		name <- "3HIV";
# 	}
# 	if(name == "3Mpa"){
# 		name <- "3MPA";
# 	}

# 	dat_temp <- tibble(
# 											pheno = dat %>% pull(pheno),
# 											log = dat %>% pull(log),
# 											int = dat %>% pull(int)
# 											);

# 	dat_temp %>% gather() -> dat_temp;

# 	i <- i + 1;
# 	p[[i]] <- eval(
# 						 dat_temp %>% 
# 							 ggplot(aes(x=value)) +
# 								geom_density(fill="black") +
# 								geom_vline(aes(xintercept=mean(pheno))) +
# 								facet_grid(. ~ key) +
# 								xlab(name) +
# 								theme_classic() +
# 								theme(
# 									axis.title.y=element_blank()
# 								)
# 						);

# }

# png(filename="_static/dens.png", width=750, height=990);
# p %>% multiplot(cols = 3);
# dev.off();

# # }}}

# }}}
# tree plot {{{

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


