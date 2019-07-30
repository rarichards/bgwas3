library("XGR")

# load data
RData.location <- "http://galahad.well.ox.ac.uk/bigdata";
xRDataLoader(RData.customised='JKscience_TS1A', RData.location=RData.location) -> dat;
# > dat %>% head()
#   ArrayAddress GeneID Symbol Naive_average LPS2_average LPS24_average ...
# 1        10008    473   RERE          9.67         8.50          9.35 ...
# 2        10037  51533   PHF7          7.31         7.19          7.21 ...
# 3        10044  55973 BCAP29          8.73         8.77          8.96 ...
# 4        10050 167153  PAPD4         10.53        10.99         10.80 ...
# 5        10068  10993    SDS          7.14         5.62          5.41 ...
# 6        10075  79724 ZNF768          7.53         7.14          7.53 ...
#
# dat %>% names()
# [1] "ArrayAddress"      "GeneID"            "Symbol"
# [4] "Naive_average"     "LPS2_average"      "LPS24_average"
# [7] "IFN_average"       "logFC_LPS2_Naive"  "fdr_LPS2_Naive"
# [10] "logFC_LPS24_Naive" "fdr_LPS24_Naive"   "logFC_INF24_Naive"
# [13] "fdr_INF24_Naive"   "logFC_LPS24_LPS2"  "fdr_LPS24_LPS2"


# create list of background genes
dat$Symbol -> background;
# > background %>% head()
# [1] "RERE"   "PHF7"   "BCAP29" "PAPD4"  "SDS"    "ZNF768"


# create list of genes significantly induced by IFN24
dat %>% 
	filter(logFC_INF24_Naive < 0 & fdr_INF24_Naive<0.01) %>%
	pull("Symbol") -> ifn24;
# > ifn24 %>% head()
# [1] "RERE"     "PAPD4"    "F3"       "LIN52"    "CD558651" "C19orf24"

xEnricherGenes(data=data, background=background, ontology="DO", ontology.algorithm="none", RData.location=RData.location) -> test;








