meso <- read.csv('meso.csv', sep=',', header=FALSE)
thermo <- read.csv('thermo.csv', sep=',', header=FALSE)
hyper <- read.csv('hyper.csv', sep=',', header=FALSE)

colnames(meso) <- c('real', 'perm')
colnames(thermo) <- c('real', 'perm')
colnames(hyper) <- c('real', 'perm')

meso_t <- t.test(meso$real, meso$perm, paired=TRUE)
thermo_t <- t.test(thermo$real, thermo$perm, paired=TRUE)
hyper_t <- t.test(hyper$real, hyper$perm, paired=TRUE)
