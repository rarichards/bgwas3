library(argparse);
#library(stringr);
#library(dplyr);
#library(readr);

parser <- ArgumentParser(description="Count unique patterns and calculate p value threshold using bonferoni");
parser$add_argument('file_patterns', help='patterns file');
#parser$add_argument('out_file', help='out file (stats)');
parser$add_argument('--alpha', default=0.05, help='family-wise error rate');

args <- parser$parse_args()
#args <- parser$parse_args(c("associations/penicillin_patterns.txt", "associations/penicillin_assoc.txt"));

statement <- paste0("LC_ALL=C sort -u -S 2014M ", args$file_patterns, " | wc -l");
unique_patterns <- strtoi(system(statement, intern=TRUE));
bonf_thresh <- args$alpha/unique_patterns;

paste0("unique_patterns\tbonf_thresh\n", unique_patterns, "\t", bonf_thresh);
# data_frame(unique_patterns=unique_patterns, bonf_thresh=bonf_thresh) %>% write_tsv(args$out_file);
