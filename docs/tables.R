library(dplyr);
library(readr);
library(stringr);
library(tidyr);

# met {{{
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

"_dat/summary.tsv" %>%
	read_tsv() %>%
	select(-bonf_thresh) %>%
	filter(significant_kmers != 0) %>%
	arrange(desc(significant_genes)) %>%
	rename(
				 "Phenotype" = pheno,
				 "Significant Kmers" = significant_kmers,
				 "Genes"= significant_genes
	) -> dat_genes;


dat_genes$Phenotype[dat_genes$Phenotype == "isoleucine_int_log"] <- "log(Isoleucine)";
dat_genes$Phenotype[dat_genes$Phenotype == "hcn_log"] <- "log(Hydrogen Cyanide)";
dat_genes$Phenotype[dat_genes$Phenotype == "ci_log"] <- "log(Ciprofloxacin)";
dat_genes$Phenotype[dat_genes$Phenotype == "at"] <- "Aztreonam";
dat_genes$Phenotype[dat_genes$Phenotype == "3HIV"] <- "3-Hydroxyisovalerate";
dat_genes$Phenotype[dat_genes$Phenotype == "3HIV_log"] <- "log(3-Hydroxyisovalerate)";
dat_genes$Phenotype[dat_genes$Phenotype == "anthranilate"] <- "Anthranilate";
dat_genes$Phenotype[dat_genes$Phenotype == "anthranilate_log"] <- "log(Anthranilate)";
dat_genes$Phenotype[dat_genes$Phenotype == "at_log"] <- "log(Aztreonam)";
dat_genes$Phenotype[dat_genes$Phenotype == "ci"] <- "Ciprofloxacin";
dat_genes$Phenotype[dat_genes$Phenotype == "co_log"] <- "log(Colistin)";
dat_genes$Phenotype[dat_genes$Phenotype == "isoleucine"] <- "Isoleucine";
dat_genes$Phenotype[dat_genes$Phenotype == "isoleucine_log"] <- "log(Isoleucine)";
dat_genes$Phenotype[dat_genes$Phenotype == "leucine"] <- "Leucine";
dat_genes$Phenotype[dat_genes$Phenotype == "leucine_log"] <- "log(Leucine)";
dat_genes$Phenotype[dat_genes$Phenotype == "methanol"] <- "Methanol";
dat_genes$Phenotype[dat_genes$Phenotype == "methanol_log"] <- "log(Methanol)";
dat_genes$Phenotype[dat_genes$Phenotype == "swarm"] <- "Swarm";
dat_genes$Phenotype[dat_genes$Phenotype == "swim_log"] <- "log(Swim)";
dat_genes$Phenotype[dat_genes$Phenotype == "tm"] <- "Tobromycin";
dat_genes$Phenotype[dat_genes$Phenotype == "tryptophan"] <- "Tryptophan";
dat_genes$Phenotype[dat_genes$Phenotype == "tryptophan_log"] <- "log(Tryptophan)";
dat_genes$Phenotype[dat_genes$Phenotype == "uracil"] <- "Uracil";
dat_genes$Phenotype[dat_genes$Phenotype == "uracil_log"] <- "log(Uracil)";

