library(phytools)
library(ape)
library(phangorn)
library(nlme)
library(phylolm)
library(ggplot2)
library(ggrepel)
library(gridExtra)
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

df$shiftedMT = 0
for(clade in unique(df$clade)) {
  this.mean = mean(df$countsMT[df$clade==clade])
  df$shiftedMT[df$clade==clade] = df$countsMT[df$clade==clade]-this.mean
}
results.df = data.frame()

# there's a problem with trait 35. we'll artificially omit it here
column.set = 1:ncol(df)
not.interesting = c(1, 2, 3, 15, 35, 56, 66, 67, 68)
cor.plot = cor.titles = list()
counter = 1
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
      mydf2 = data.frame(label=df$Scientific.Name, x=ifelse(df[,positive.column]==positive.label, 1, 0), y=df$shiftedMT, clade=df$clade)
        # loop through leaves in taxonomy tree
      rownames(mydf2) = mydf2$label
      mydf2 = mydf2[which(mydf2$label %in% positive.clade$tip.label),]
        # construct data items for use below
        #nametips = positive.clade$tip.label
        #fx = as.factor(x)			
        #mydf = data.frame(y,x)			
        #mydf2 = data.frame(label=nametips, x=fx, y=y)

      expt.label = paste(c(colnames(df)[positive.column], positive.label), collapse=" ")
      
    
      
#      cor.test(mydf2$x, mydf2$y)
        #### PLM
          print(paste(c(
            "Trying", positive.column, colnames(df)[positive.column], positive.label, 
            "for", length(positive.tips), "positives and positive clade of size", 
            length(positive.clade$tip.label)), collapse=" "))
        if(all(mydf2$x=="0") | all(mydf2$x=="1") | all(mydf2$y==mydf2$y[1])) {
          coef = 0
          pval = 1
          plm.coef = 0
          plm.pval = 1
        } else {
          plm = phylolm(y~x, mydf2, positive.clade, model="BM")
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
          
          toplot = c()
          for(clade in unique(mydf2$clade)) {
            if(length(unique(mydf2$x[mydf2$clade==clade])) > 1) { toplot = c(toplot, clade) }
          }
          if(is.na(plm.pval)) {
            plot.label = paste(c(expt.label, ", p=NAN"), collapse="")
          } else  if(plm.pval == 0) {
            plot.label = paste(c(expt.label, ", p<1e-50"), collapse="")
          } else {
            plot.label = paste(c(expt.label, ", p=", signif(plm.pval, digits=2)), collapse="")
          }
          toplot.df = mydf2[mydf2$clade %in% toplot,]
          cor.plot[[counter]] = ggplot(toplot.df, aes(x=clade,y=y,color=factor(x))) + # geom_boxplot(position="dodge") + 
            geom_quasirandom(alpha=1,dodge.width=0.5) +
            theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
            ggtitle(plot.label) + ylab("Count difference") + labs(color="Presence")
          counter =counter+1
          
          print(paste(c("Done ", positive.column, positive.label, length(positives))))
          results.df = rbind(results.df, data.frame(col=positive.column, colname=colnames(df)[positive.column],
                                                    positive.label=positive.label,
                                                    plm.coef=plm.coef,plm.pval=plm.pval,
                                                    n.positive=length(positives), clade.positive=length(positive.clade$tip.label),
                                                    mrca.positive=tree.labels[MRCA], mincount=mincount))
        write.csv(results.df, "test-clademean.csv", quote=FALSE, row.names=FALSE)
        #write(paste(positive.column, positive.label, coef, pval, length(positives), length(positive.clade$tip.label)), file="23a-mt-tmp.txt", sep=",", append=TRUE)
         
      }
    }
  }
}

valid.df = results.df[results.df$mincount > 3 & results.df$n.positive > 6 & (results.df$plm.pval < 0.05) ,]

p.label = function(p1, n) {
  if(p1 < 0.05/n ) { return("**")}
 
  if(p1 < 0.05) { return("*")}
  return("-/-")
}

valid.df$label = ""
valid.df$p.cat = 0
for(i in 1:nrow(valid.df)) {
  
  valid.df$label[i] = paste(c(valid.df$colname[i], ":\n", valid.df$positive.label[i]), collapse="")
  valid.df$p.cat[i] = p.label(valid.df$plm.pval[i], nrow(results.df))
}
valid.df$p.cat = factor(valid.df$p.cat, levels=c("**", "*", "-"))

g.mt.plm.clademean = ggplot(valid.df, aes(x=plm.coef, y=log(-log(plm.pval)), label=label, color=p.cat)) + 
  geom_point() + geom_text_repel(max.overlaps=50, size=3, lineheight=0.75) + 
  theme_light() + labs(title="Shifted MT gene count PLM", x ="PLM coefficient", y = "log(-log(p))", color="p profile")

g.mt.plm.clademean

sf = 2
png("mt-test-clademean.png", width=400*sf, height=300*sf, res=72*sf)
g.mt.plm.clademean
dev.off()

which(results.df$colname == "locomotion")
which(results.df$colname == "woodiness")
which(results.df$colname == "plant.growth.form")
which(results.df$colname == "parasiteOf")
which(results.df$colname == "desert")

cor.plot[[48]]

png("mt-test-clademean-specifics.png", width=600, height=400)
grid.arrange(cor.plot[[12]], cor.plot[[122]], cor.plot[[22]], cor.plot[[8]], cor.plot[[165]], nrow=2)
dev.off()
