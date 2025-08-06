library(ape)
library(nlme)
library(geiger)
library(phangorn)

data <- read.csv('metadata.csv', header = TRUE, row.names = 1)
tree <- read.nexus('tree.nex')
obj <- name.check(tree, data)
pruned_tree <- drop.tip(tree, obj$tree_not_data)
mid_tree <- midpoint(pruned_tree)
data$Genus <- rownames(data)
data <- data[data$Genus %in% mid_tree$tip.label, ]
bm <- corBrownian(phy = mid_tree, form = ~Genus)

model1 <- gls(Growth_Temp ~ Genome_GC, data = data, correlation = bm)
model2 <- gls(Growth_Temp ~ CDS_room, data = data, correlation = bm)
model3 <- gls(Growth_Temp ~ Core_room, data = data, correlation = bm)
