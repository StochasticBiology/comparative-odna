library(phytools)
library(ape)
library(phangorn)
library(nlme)
library(phylolm)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(rcompanion)
library(ggbeeswarm)

# this function converts a species name string from the Newick format which Common Taxonomy Tree gives us into a simpler lower-case, no quotes version comparable to Kostas' dataset
convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

# read phylogeny previously downloaded from Common Taxonomy Tree. this has previously been cleaned:
# (i) all one line; (ii) extra set of brackets wrap the whole entity; (iii) nodes names don't contain special characters
tree = read.newick("tree-for-traits-clean-mt.phy")
#my.correlation = corBrownian(phy=tree)
tree$tip.label = convname(tree$tip.label)
tree.labels = c(tree$tip.label, tree$node.label)
root = which(tree.labels=="Eukaryota")
clade.refs = Children(tree, root)
clade.names = tree.labels[clade.refs]

# read Kostas' dataset
df = read.table("MTFull22.txt", sep="\t", header=T, stringsAsFactors = TRUE)

# manually fix bugs
df$plant.growth.form[which(df$Scientific.Name=="triodia sylvina")] = NA
df$life.cycle.habit[grep("lucilia", df$Scientific.Name)] = NA

# the dataset is read as factors by default -- we'll awkwardly convert between data types
df$Scientific.Name = as.character(df$Scientific.Name)
df$Ancestry = as.character(df$Ancestry)

df$clade = ""
to.drop = c()
for(i in 1:nrow(df)) {
  tip.ref = which(tree.labels == df$Scientific.Name[i])
  if(length(tip.ref) > 0) {
    ancestors = Ancestors(tree, tip.ref)
    ancestor = ancestors[which(ancestors %in% clade.refs)]
    df$clade[i] = tree.labels[ancestor]
  } else {
    to.drop = c(to.drop, i)
  }
}

df = df[-to.drop,]
#df = df[!is.na(df$y),]

#df = df[df$clade != "",]
results.df = data.frame()

# there's a problem with trait 35. we'll artificially omit it here
# we're omitting more here, including factors with many many levels
column.set = 1:ncol(df)
not.interesting = c(1, 2, 3, 14, 15, 16, 28, 31, 35, 43, 56, 63, 64, 65, 66, 67, 68, 76, 78)
cor.plot = cor.titles = list()
counter = 1
for(positive.column in column.set[-not.interesting]) {
  expt.label = colnames(df)[positive.column]
  print(paste(c(positive.column, expt.label), collapse=" "))
  mydf2 = data.frame(label=df$Scientific.Name, x=df[,positive.column], y=df$countsMT, clade=df$clade)
  rownames(mydf2) = mydf2$label
  
  mydf2$x[is.na(mydf2$x)] = ""
  if(length(unique(mydf2$x)) > 1) {
    if(length(unique(mydf2$x)) > 2) {
      mydf2 = mydf2[mydf2$x != "" & !is.na(mydf2$x),]
    }
    if(length(unique(mydf2$clade)) > 1) {
      test.SRH = scheirerRayHare(y ~ x + clade, data=mydf2)
    } else {
      test.SRH = kruskal.test(mydf2$y, mydf2$x)
    }
    srh.pval = test.SRH$p.value[1]
  } else {
    srh.pval = 1
  }
  
  if(is.na(srh.pval)) {
    plot.label = paste(c(expt.label, ", p=NAN"), collapse="")
  } else  if(srh.pval == 0) {
    plot.label = paste(c(expt.label, ", p<1e-50"), collapse="")
  } else {
    plot.label = paste(c(expt.label, ", p=", signif(srh.pval, digits=2)), collapse="")
  }
  
  toplot = c()
  for(clade in unique(mydf2$clade)) {
    if(length(unique(mydf2$x[mydf2$clade==clade])) > 1) { toplot = c(toplot, clade) }
  }
  toplot.df = mydf2[mydf2$clade %in% toplot,]
  cor.plot[[counter]] = ggplot(toplot.df, aes(x=clade,y=y,color=x)) + # geom_boxplot(position="dodge") + 
    geom_quasirandom(alpha=1,dodge.width=0.5) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    ggtitle(plot.label) 
  counter=counter+1
  results.df = rbind(results.df, data.frame(label=expt.label,col=positive.column, srh.pval=srh.pval))
}

valid.df = results.df[!is.na(results.df$srh.pval) & results.df$srh.pval < 0.05/nrow(results.df),]

# for example
cor.plot[[7]]
to.plot = cor.plot[as.numeric(row.names(valid.df))]
sf=2
png("mt-blocking.png", width=1000*sf, height=1000*sf, res=72*sf)
grid.arrange(grobs=to.plot)
dev.off()
