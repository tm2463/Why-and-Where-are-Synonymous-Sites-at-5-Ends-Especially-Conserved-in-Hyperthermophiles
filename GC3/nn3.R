library(ggplot2)
library(tidyverse)

nn3_bias <- read.csv('nn3.csv', sep=',')
rownames(nn3_bias) <- nn3_bias$X
nn3_bias$X <- NULL

nn3 <- nn3_bias %>%
	rownames_to_column(var = "Group") %>%
	pivot_longer(cols = -Group, names_to = "nt", values_to = "rscu")

heatmap <- ggplot(nn3, aes(x=nt, y=Group, fill=rscu)) +
	geom_tile() +
	coord_fixed() +
	scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 1) +
	labs(title='NN3 Bias', x=NULL, y=NULL)

ggsave('nn3.png', plot=heatmap, width=4, height=4)