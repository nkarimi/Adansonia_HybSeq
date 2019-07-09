library(reshape2)
library(tidyr)
library(lattice)
library(dplyr)
library(ggplot2)



count.table <- read.table("count.table", sep="\t", header=F)
count.uniq <- unique(count.table)

count.table.rs=as.data.frame(acast(count.uniq, V1~V2, value.var="V3"))
count.table.rs=data.matrix(count.table.rs, rownames.force=T)
count.table.rs[is.na(count.table.rs)] <- 0 

count.table.pp=sweep(count.table.rs[c(1:3),], 1, 4, FUN="/")
count.table.di=sweep(count.table.rs[c(4:nrow(count.table.rs)),], 1, 2, FUN="/")
m=rbind(count.table.pp, count.table.di)

md <- melt(m)

png("alleles.png", 9000, 3000, pointsize=12, res=600)
ggplot(md, aes(x = Var2, y = Var1, fill = factor(value))) + geom_tile() + scale_fill_manual(values=c("red", "white", "green3", "blue", "black"),labels=c("0 alleles","1 allele","2 alleles","3 alleles","4 alleles"),name="Number of alleles\nrecovered") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("Taxon")
dev.off()



#white = 0 recovered
#black = 2alleles (diploids) and 4alleles (tetraploid)
#green = 1 allele (diploids) and 2 allele (tetraploids)
#light blue = 1 allele for tetraploid
#blue = 3 alleles (tetraploid)
