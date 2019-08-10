library(argparse);
library(dplyr);
library(readr);

parser <- ArgumentParser(description="count genes");
parser$add_argument('pvalues', help="list of p values");
parser$add_argument('--prefix', help="", default="plot");
parser$add_argument('--confidence-interval', help="", default=0.95);

args <- parser$parse_args();
