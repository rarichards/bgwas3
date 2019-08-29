suppressWarnings(suppressMessages(library(argparse)));
suppressWarnings(suppressMessages(library(dplyr)));
suppressWarnings(suppressMessages(library(readr)));
suppressWarnings(suppressMessages(library(ggplot2)));

get_args <- function(){
	description <- "plot a density plot of a continuous variable";
	parser <- ArgumentParser(description = description); 
	parser$add_argument("input", help='tsv file where a column is a variable to be plotted');
	parser$add_argument('output', help='path of png file output');
	parser$add_argument('--column', help='column to use', type="integer", default=1);
	parser$add_argument('--width', help='width (px)', type="integer", default=300);
	parser$add_argument('--height', help='height (px)', type="integer", default=300);
	parser$add_argument('--axis', help='axis (0=false, 1=true)', type="integer", default=1);
	parser$add_argument('--xaxis', help='x axis title', default=NULL);

	return(parser$parse_args());
}

main <- function(dat, png_path, xaxis, width, height, axis){
	if(axis){
		if(is.null(xaxis)){
			xaxis <- names(dat)[1];
		}

		dat %>%
			ggplot(aes_string(x=names(dat)[1])) +
				geom_density(fill="black") +
				ylab("Density") +
				xlab(xaxis) +
				theme_classic() -> p;
	}else{
		dat %>%
			ggplot(aes_string(x=names(dat)[1])) +
				geom_density(fill="black") +
				theme_classic() +
				theme(
					axis.title=element_blank(),
					axis.text=element_blank(),
					axis.line=element_blank(),
					axis.ticks=element_blank()
				) -> p;
	}

	# dpi = 300;
	# width = width/dpi;
	# height = height/dpi

	ggsave(plot=p, filename=png_path, width=width, height=height);
	# ggsave(plot=p, filename=png_path, width=width, height=height, units="in", dpi=dpi);
}

args <- get_args();

args$input %>% 
	read_tsv %>%
	select(args$column) -> dat;

axis = args$axis == 1;

main(dat, args$output, args$xaxis, args$width, args$height, axis);
