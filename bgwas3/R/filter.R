library(argparse);
library(stringr);
library(dplyr);

parser <- ArgumentParser(description="Count unique patterns, filter associated kmers with bonferoni and generate a stats file");
parser$add_argument('file_patterns', help='patterns file');
parser$add_argument('file_assoc', help='associations file');
parser$add_argument('--alpha', default=0.05, help='Family-wise error rate');

args <- parser$parse_args()
#args <- parser$parse_args(c("associations/penicillin_patterns.txt", "associations/penicillin_assoc.txt"));

statement <- paste0("LC_ALL=C sort -u -S 2014M ", args$file_patterns, " | wc -l");
unique_patterns <- strtoi(system(statement, intern=TRUE));
bonf_thresh <- args$alpha/unique_patterns;

file_filtered <- args$file_assoc %>% str_sub(1, -5) %>% paste0("_filtered.txt");

statement <- paste0("cat <(head -1 ",  args$file_assoc,  ") <(awk '$4<", bonf_thresh, "{print $0}' ", args$file_assoc, ") > ", file_filtered);
system(statement);
