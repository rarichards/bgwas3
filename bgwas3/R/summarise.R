library(argparse);
library(dplyr);
library(readr);

parser <- ArgumentParser(description="Count unique patterns and calculate p value threshold using bonferoni");
parser$add_argument('file_patterns', help='patterns file');
parser$add_argument('out_file', help='out file (stats)');
parser$add_argument('--alpha', default=0.05, help='family-wise error rate');

args <- parser$parse_args()

