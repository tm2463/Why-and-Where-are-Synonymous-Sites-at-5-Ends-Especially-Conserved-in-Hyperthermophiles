library(ape)
library(nlme)
library(geiger)
library(phangorn)
library(ggplot2)

data <- read.csv('metadata2.csv', header = TRUE, row.names = 1)
tree <- read.nexus('tree.nex')
obj <- name.check(tree, data)
pruned_tree <- drop.tip(tree, obj$tree_not_data)
mid_tree <- midpoint(pruned_tree)
data$Genus <- rownames(data)
data <- data[data$Genus %in% mid_tree$tip.label, ]
bm <- corBrownian(phy = mid_tree, form = ~Genus)

model1 <- gls(Growth_Temp ~ forty, data = data, correlation = bm)
table <- summary(model1)$tTable
p1 <- table["forty", "p-value"]

model2 <- gls(Growth_Temp ~ fifty, data = data, correlation = bm)
table <- summary(model2)$tTable
p2 <- table["fifty", "p-value"]

model3 <- gls(Growth_Temp ~ sixty, data = data, correlation = bm)
table <- summary(model3)$tTable
p3 <- table["sixty", "p-value"]

model4 <- gls(Growth_Temp ~ seventy, data = data, correlation = bm)
table <- summary(model4)$tTable
p4 <- table["seventy", "p-value"]

model5 <- gls(Growth_Temp ~ eighty, data = data, correlation = bm)
table <- summary(model5)$tTable
p5 <- table["eighty", "p-value"]

model6 <- gls(Growth_Temp ~ ninety, data = data, correlation = bm)
table <- summary(model6)$tTable
p6 <- table["ninety", "p-value"]

p_values <- c(forty = p1, fifty = p2, sixty = p3, seventy = p4, eighty = p5, ninety = p6)
temps <- c(40, 50, 60, 70, 80, 90)

df <- data.frame(Temperature=temps, p_value=as.numeric(p_values))

line <- ggplot(data=df, aes(x=Temperature, y=p_value)) +
	geom_line() +
	geom_point() +
	geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  	labs(y = "p-value", title = "Temp vs PGLS p-value")

ggsave('pgls.png', plot=line, width=6, height=4)
