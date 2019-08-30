library(dplyr);
library(readr);
library(stringr);
library(tidyr);
library(ggplot2);

"phenos.tsv" %>% 
	read_tsv() %>%
	select(starts_with("pheno")) -> dat;

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
