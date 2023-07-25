library(phytools)
library(ape)
library(phangorn)
library(nlme)
library(phylolm)
library(ggplot2)
library(ggrepel)
library(gridExtra)

# this function converts a species name string from the Newick format which Common Taxonomy Tree gives us into a simpler lower-case, no quotes version comparable to Kostas' dataset
convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

# read phylogeny previously downloaded from Common Taxonomy Tree. this has previously been cleaned:
# (i) all one line; (ii) extra set of brackets wrap the whole entity; (iii) nodes names don't contain special characters
tree = read.newick("tree-for-traits-clean-pt.phy")
#my.correlation = corBrownian(phy=tree)
tree$tip.label = convname(tree$tip.label)
tree.labels = c(tree$tip.label, tree$node.label)
# assign clade labels
root = which(tree.labels=="Eukaryota")
clade.refs = Children(tree, root)
clade.names = tree.labels[clade.refs]

# read Kostas' dataset
df = read.table("PTFull22.txt", sep="\t", header=T, stringsAsFactors = TRUE)

# manually fix bugs
df$plant.growth.form[which(df$Scientific.Name=="triodia sylvina")] = NA
df$life.cycle.habit[grep("lucilia", df$Scientific.Name)] = NA

# the dataset is read as factors by default -- we'll awkwardly convert between data types
df$Scientific.Name = as.character(df$Scientific.Name)
df$Ancestry = as.character(df$Ancestry)

# either assign clade label to dataframe rows, or drop if not found in phylogeny
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

results.df = data.frame()

# there's a problem with trait 35. we'll artificially omit it here
column.set = 1:ncol(df)
not.interesting = c(1, 2, 3, 15, 35, 56, 66, 67, 68)
cor.plot = cor.titles = list()
for(positive.column in column.set[-not.interesting]) {
  label.set = unique(df[,positive.column])
  for(positive.label in label.set) {
    if(is.factor(label.set) & positive.label != "" & !is.na(positive.label)) {
      # get set of positives names
      positives = df[df[,positive.column]==positive.label,1]
      positive.tips = which(tree$tip.label %in% positives)
      if(length(positive.tips) > 1) {
        MRCA = findMRCA(tree, positive.tips)
        positive.clade = extract.clade(tree, MRCA)
        pc.labels = c(positive.clade$tip.label, positive.clade$node.label)
        pc.root = which(pc.labels == tree.labels[MRCA])
        positive.clade = root(positive.clade, node=pc.root)
        # now construct a dataframe for use in other approaches
      mydf2 = data.frame(label=df$Scientific.Name, x=ifelse(df[,positive.column]==positive.label, 1, 0), y=df$countsPT)
        # loop through leaves in taxonomy tree
      rownames(mydf2) = mydf2$label
      mydf2 = mydf2[which(mydf2$label %in% positive.clade$tip.label),]
        # construct data items for use below
        #nametips = positive.clade$tip.label
        #fx = as.factor(x)			
        #mydf = data.frame(y,x)			
        #mydf2 = data.frame(label=nametips, x=fx, y=y)

      expt.label = paste(c(colnames(df)[positive.column], positive.label), collapse=" ")
      cor.plot[[length(cor.plot)+1]] = ggplot(mydf2, aes(x=x,y=y)) + geom_jitter(width=0.2) + 
        geom_smooth(method="lm") + ggtitle(expt.label) + theme(plot.title = element_text(size = 12))
      cor.titles[[length(cor.titles)+1]] = expt.label
      
#      cor.test(mydf2$x, mydf2$y)
        #### PGLM
          print(paste(c(
            "Trying", positive.column, colnames(df)[positive.column], positive.label, 
            "for", length(positive.tips), "positives and positive clade of size", 
            length(positive.clade$tip.label)), collapse=" "))
        if(all(mydf2$x=="0") | all(mydf2$x=="1") | all(mydf2$y==mydf2$y[1])) {
          coef = 0
          pval = 1
          pglm.coef = 0
          pglm.pval = 1
          plm.coef = 0
          plm.pval = 1
        } else {
          pglm = phyloglm(y~x, mydf2, positive.clade, method="poisson_GEE")
          plm = phylolm(y~x, mydf2, positive.clade, model="BM")
          pglm.coef = summary(pglm)$coefficients[2,1]
          pglm.pval = summary(pglm)$coefficients[2,4]
          plm.coef = summary(plm)$coefficients[2,1]
          plm.pval = summary(plm)$coefficients[2,4]
        }
         
          
          mincount = Inf
          if(length(unique(mydf2$y[mydf2$x==1])) < 3 | length(unique(mydf2$y[mydf2$x==0])) < 3) {
          for(i in unique(mydf2$y)) {
            c = length(which(mydf2$y[mydf2$x==1] == i))
            if(c < mincount) { mincount = c }
            c = length(which(mydf2$y[mydf2$x==0] == i))
            if(c < mincount) { mincount = c }
          }
          }
          
          print(paste(c("Done ", positive.column, positive.label, length(positives))))
          results.df = rbind(results.df, data.frame(col=positive.column, colname=colnames(df)[positive.column],
                                                    positive.label=positive.label,pglm.coef=pglm.coef,pglm.pval=pglm.pval,
                                                    plm.coef=plm.coef,plm.pval=plm.pval,
                                                    n.positive=length(positives), clade.positive=length(positive.clade$tip.label),
                                                    mrca.positive=tree.labels[MRCA], mincount=mincount))
        write.csv(results.df, "test.csv", quote=FALSE, row.names=FALSE)
        #write(paste(positive.column, positive.label, coef, pval, length(positives), length(positive.clade$tip.label)), file="23a-pt-tmp.txt", sep=",", append=TRUE)
         
      }
    }
  }
}



valid.df = results.df[results.df$mincount > 3 & (results.df$pglm.pval < 0.05 | results.df$plm.pval < 0.05) ,]
valid.df = results.df[results.df$mincount > 3 & results.df$n.positive > 6 & (results.df$pglm.pval < 0.05 | results.df$plm.pval < 0.05) ,]

p.label = function(p1, p2, n) {
  if(p1 < 0.05/n & p2 < 0.05/n) { return("**/**")}
  if((p1 < 0.05/n & p2 < 0.05) | (p1 < 0.05 & p2 < 0.05/n)) { return("**/*")}
  if((p1 < 0.05/n | p2 < 0.05/n)) { return("**/-")}
  if((p1 < 0.05 & p2 < 0.05)) { return("*/*")}
  if((p1 < 0.05 | p2 < 0.05)) { return("*/-")}
  return("-/-")
}

valid.df$label = ""
valid.df$p.cat = 0
for(i in 1:nrow(valid.df)) {
  
  valid.df$label[i] = paste(c(valid.df$colname[i], ":\n", valid.df$positive.label[i], "\n", valid.df$mrca.positive[i]), collapse="")
  valid.df$p.cat[i] = p.label(valid.df$plm.pval[i], valid.df$pglm.pval[i], nrow(results.df))
}
valid.df$p.cat = factor(valid.df$p.cat, levels=c("**/**", "**/*", "**/-", "*/*", "*/-", "-/-"))
g.pglm = ggplot(valid.df, aes(x=pglm.coef, y=log(-log(pglm.pval)), label=label, color=p.cat)) + 
  geom_point() + geom_text_repel(max.overlaps=50, size=2) + xlim(NA,0.8)
g.plm = ggplot(valid.df, aes(x=plm.coef, y=log(-log(plm.pval)), label=label, color=p.cat)) + 
  geom_point() + geom_text_repel(max.overlaps=50, size=2) + xlim(NA,40)
grid.arrange(g.pglm, g.plm, nrow=2)

sf = 2
png("pt-test.png", width=600*sf, height=400*sf, res=72*sf)
grid.arrange(g.pglm, g.plm)
dev.off()

#results.df[which(results.df$colname == "trophic.guild.x"),]
#cor.plot[[9]]
#cor.plot[[11]]
#results.df[which(results.df$colname == "unicellularity"),]
#cor.plot[[118]]
#results.df[which(results.df$colname == "ecomorphological.guild"),]
#cor.plot[[45]]
#results.df[which(results.df$colname == "parasiteOf"),]
#cor.plot[[6]]

png("pt-test-specifics.png", width=600, height=400)
grid.arrange(cor.plot[[9]], cor.plot[[11]], cor.plot[[118]], cor.plot[[45]], cor.plot[[6]], nrow=2)
dev.off()
