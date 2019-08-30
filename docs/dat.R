library(dplyr);
library(readr);
library(stringr);
library(tidyr);

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

