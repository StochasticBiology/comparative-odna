library(phytools)
library(ape)
library(phangorn)
library(nlme)
library(phylolm)

# this function converts a species name string from the Newick format which Common Taxonomy Tree gives us into a simpler lower-case, no quotes version comparable to Kostas' dataset
convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

# read phylogeny previously downloaded from Common Taxonomy Tree. this has previously been cleaned:
# (i) all one line; (ii) extra set of brackets wrap the whole entity; (iii) nodes names don't contain special characters
tree = read.newick("tree-for-traits-clean-mt.phy")
tree.tips = convname(tree$tip.label)
tree.labels = c(tree$tip.label, tree$node.label)

metazoa.label = which(tree.labels=="Metazoa")
metazoa.tips = Descendants(tree, metazoa.label, type="tips")[[1]]
tree.no.metazoa = drop.tip(tree, metazoa.tips)
tnm.labels = c(tree.no.metazoa$tip.label, tree.no.metazoa$node.label)
tnm.root = which(tnm.labels == "Eukaryota")
tree.no.metazoa = root(tree.no.metazoa, node=tnm.root)
tree = tree.no.metazoa

#my.correlation = corBrownian(phy=tree)
tree$tip.label = convname(tree$tip.label)
tree.labels = c(tree$tip.label, tree$node.label)

# read Kostas' dataset
df = read.table("MTFull22.txt", sep="\t", header=T, stringsAsFactors = TRUE)

# the dataset is read as factors by default -- we'll awkwardly convert between data types
df$Scientific.Name = as.character(df$Scientific.Name)
df$Ancestry = as.character(df$Ancestry)

results.df = data.frame()

# there's a problem with trait 35. we'll artificially omit it here
column.set = 1:ncol(df)
not.interesting = c(1, 2, 3, 15, 35, 56, 66, 67, 68)
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
        # now construct a dataframe for use in other approaches
      mydf2 = data.frame(label=df$Scientific.Name, x=ifelse(df[,positive.column]==positive.label, 1, 0), y=df$countsMT)
        # loop through leaves in taxonomy tree
      rownames(mydf2) = mydf2$label
      mydf2 = mydf2[which(mydf2$label %in% positive.clade$tip.label),]
        # construct data items for use below
        #nametips = positive.clade$tip.label
        #fx = as.factor(x)			
        #mydf = data.frame(y,x)			
        #mydf2 = data.frame(label=nametips, x=fx, y=y)

        #### PGLM
          print(paste(c(
            "Trying", positive.column, colnames(df)[positive.column], positive.label, 
            "for", length(positive.tips), "positives and positive clade of size", 
            length(positive.clade$tip.label)), collapse=" "))
        if(all(mydf2$x=="0") | all(mydf2$x=="1") | all(mydf2$y==mydf2$y[1])) {
          coef = 0
          pval = 1
        } else {
          pglm = phyloglm(y~x, mydf2, positive.clade, method="poisson_GEE")
          coef = summary(pglm)$coefficients[2,1]
          pval = summary(pglm)$coefficients[2,4]
       
        }
          print(paste(c("Done ", positive.column, positive.label, coef, pval, length(positives))))
          results.df = rbind(results.df, data.frame(col=positive.column, colname=colnames(df)[positive.column],
                                                    positive.label=positive.label,coef=coef,pval=pval,
                                                    n.positive=length(positives), clade.positive=length(positive.clade$tip.label),
                                                    mrca.positive=tree.labels[MRCA]))
        write.csv(results.df, "mac-results-df-mt-no-metazoa-pglm.csv", quote=FALSE, row.names=FALSE)
        #write(paste(positive.column, positive.label, coef, pval, length(positives), length(positive.clade$tip.label)), file="23a-mt-tmp.txt", sep=",", append=TRUE)
         
      }
    }
  }
}

