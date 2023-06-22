library(phytools)
library(phyreg)
library(ape)
library(phangorn)
library(geiger)
library(nlme)
library(caper)

# this function converts a species name string from the Newick format which Common Taxonomy Tree gives us into a simpler lower-case, no quotes version comparable to Kostas' dataset
convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

# read phylogeny previously downloaded from Common Taxonomy Tree. this has previously been cleaned:
# (i) all one line; (ii) extra set of brackets wrap the whole entity; (iii) nodes names don't contain special characters
tree = read.newick("tree-for-traits-clean-mt.phy")
#my.correlation = corBrownian(phy=tree)
tree.tips = convname(tree$tip.label)
tree.labels = c(tree$tip.label, tree$node.label)

# read Kostas' dataset
df = read.table("MTFull22.txt", sep="\t", header=T, stringsAsFactors = TRUE)

# the dataset is read as factors by default -- we'll awkwardly convert between data types
df$Scientific.Name = as.character(df$Scientific.Name)
df$Ancestry = as.character(df$Ancestry)

results.df = data.frame()

column.set = 1:ncol(df)
not.interesting = c(1, 2, 3, 15, 35, 56, 66, 67, 68)
for(positive.column in column.set[-not.interesting]) {
  label.set = unique(df[,positive.column])
  for(positive.label in label.set) {
    if(is.factor(label.set) & positive.label != "" & !is.na(positive.label)) {
      # get set of positives names
      positives = df[df[,positive.column]==positive.label,1]
      positive.tips = which(tree.tips %in% positives)
      if(length(positive.tips) > 1) {
        MRCA = findMRCA(tree, positive.tips)
        positive.clade = extract.clade(tree, MRCA)
        # now construct a dataframe for use in other approaches
      mydf2 = data.frame()
        # loop through leaves in taxonomy tree
        
          # get a name to cross-reference with Kostas' dataset
          thesenames = positive.clade$tip.label
          print("Building df")
          for(i in 1:length(thesenames)) {
            cname = convname(thesenames[i])
                ref = which(df$Scientific.Name == cname)
            if(length(ref) > 0) {
              if(cname %in% positives) { this.x = 1 } else { this.x = 0 }
              mydf2 = rbind(mydf2, data.frame(label=thesenames[i], x=as.factor(this.x), y=df$countsMT[ref]))
            }
          }
            
            
        # construct data items for use below
        #nametips = positive.clade$tip.label
        #fx = as.factor(x)			
        #mydf = data.frame(y,x)			
        #mydf2 = data.frame(label=nametips, x=fx, y=y)

        #### PGLS via Luke Harmon https://lukejharmon.github.io/ilhabela/instruction/2015/07/03/PGLS/
        # this eats memory (at least 2GB) and takes perhaps dozens of minutes
        # name.check(mydf2, tree)
          print(paste(c(
            "Trying", positive.column, colnames(df)[positive.column], positive.label, 
            "for", length(positive.tips), "positives and positive clade of size", 
            length(positive.clade$tip.label)), collapse=" "))
        if(all(mydf2$x=="0") | all(mydf2$x=="1") | all(mydf2$y==mydf2$y[1])) {
          coef = 0
          pval = 1
        } else {
        my.correlation = corBrownian(phy=positive.clade, form=~label)
        pglsModel2<-gls(y~x, correlation=my.correlation, data=mydf2, method="ML", na.action=na.omit)
        coef = coef(pglsModel2)[2]
        pval = anova(pglsModel2)$'p-value'[2]
        }
          print(paste(c("Done ", positive.column, positive.label, coef, pval, length(positives))))
          results.df = rbind(results.df, data.frame(col=positive.column, colname=colnames(df)[positive.column],
                                                    positive.label=positive.label,coef=coef,pval=pval,
                                                    n.positive=length(positives), clade.positive=length(positive.clade$tip.label),
                                                    mrca.positive=tree.labels[MRCA]))
        write.csv(results.df, "mac-results-df-mt.csv", quote=F, row.names=FALSE)
        #write(paste(positive.column, positive.label, coef, pval, length(positives), length(positive.clade$tip.label)), file="23a-mt-tmp.txt", sep=",", append=TRUE)
         
      }
    }
  }
}

