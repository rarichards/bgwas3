library("dplyr")
library("readr")
library("tidyr")
library("argparse")

parser <- ArgumentParser(description="Seperate one tsv file into multiple individual files")
parser$add_argument("file", help="tsv file")
parser$add_argument("--id", default="1", help="id column index")
args = parser$parse_args()

args$fils %>% 
		read_tsv() %>% 
		rename("id" = args$id) -> dat;

for(col in colnames(dat)){
	if(col !== "id"){
		dat %>% 
			select(c("id", col)) %>%
			write_tsv(path = col.tsv);
	}
}
