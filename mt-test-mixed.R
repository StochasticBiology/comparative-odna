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
      mydf2 = data.frame(label=df$Scientific.Name, x=ifelse(df[,positive.column]==positive.label, 1, 0), y=df$countsMT, clade=df$clade)
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
          
          # "PMIXED"
          if(length(unique(mydf2$clade)) == 1) {
            lmm.coef = glmm.coef = 0
            lmm.pval = glmm.pval = 1
            mod.lmm = lm(y ~ x, data=mydf2)
            best.lmm = summary(mod.lmm)
            if(nrow(best.lmm$coefficients) > 1) {
            lmm.coef = best.lmm$coefficients[2,1]
            lmm.pval = best.lmm$coefficients[2,4]
            }
            mod.glmm = glm(y ~ x, data=mydf2, family=poisson())
            best.glmm = summary(mod.glmm)
            if(nrow(best.glmm$coefficients) > 1) {
              glmm.coef = best.glmm$coefficients[2,1]
              glmm.pval = best.glmm$coefficients[2,4]
            }
            
          } else {
            ### SINGULARITY ISSUES HERE
            tryCatch( {
              mod.lmm = lme(y ~ x, random = ~ x | clade, data=mydf2, method="ML", control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
              mod.lmm1 = lme(y ~ x, random = ~ 1 | clade, data=mydf2, method="ML", control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
              if(AIC(mod.lmm) < AIC(mod.lmm1)) {
                best.lmm = summary(mod.lmm) 
              } else {
                best.lmm = summary(mod.lmm1)
              }
              lmm.coef = best.lmm$tTable[2,1]
              lmm.pval = best.lmm$tTable[2,5]
            },  error = function(e) {
              lmm.coef = 0
              lmm.pval = 1
            }) 
         
            tryCatch( {
          mod.glmm = glmer(y ~ x + (x | clade), data=mydf2, poisson)
          mod.glmm1 = glmer(y ~ x + (1 | clade), data=mydf2, poisson)
          if(AIC(mod.glmm) < AIC(mod.glmm1)) {
            best.glmm = summary(mod.glmm) 
          } else {
            best.glmm = summary(mod.glmm1)
          }
          glmm.coef = best.glmm$coefficients[2,1]
          glmm.pval = best.glmm$coefficients[2,4]
            },  error = function(e) {
              glmm.coef = 0
              glmm.pval = 1
            }) 
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
                                                    lmm.coef=lmm.coef,lmm.pval=lmm.pval,
                                                    glmm.coef=glmm.coef,glmm.pval=glmm.pval,
                                                    n.positive=length(positives), clade.positive=length(positive.clade$tip.label),
                                                    mrca.positive=tree.labels[MRCA], mincount=mincount))
        write.csv(results.df, "test.csv", quote=FALSE, row.names=FALSE)
        #write(paste(positive.column, positive.label, coef, pval, length(positives), length(positive.clade$tip.label)), file="23a-mt-tmp.txt", sep=",", append=TRUE)
         
      }
    }
  }
}

results.df[!is.na(results.df$lmm.pval) & results.df$lmm.pval < 0.05/nrow(results.df),]
results.df[!is.na(results.df$glmm.pval) & results.df$glmm.pval < 0.05/nrow(results.df),]

results.df[!is.na(results.df$glmm.pval) & results.df$glmm.pval < 0.05,]

valid.df = results.df[results.df$mincount > 3 & (results.df$glmm.pval < 0.05 | results.df$lmm.pval < 0.05) ,]
valid.df = results.df[results.df$mincount > 3 & results.df$n.positive > 6 & (results.df$glmm.pval < 0.05 | results.df$lmm.pval < 0.05) ,]

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
  
  valid.df$label[i] = paste(c(valid.df$colname[i], ":\n", valid.df$positive.label[i]), collapse="")
  valid.df$p.cat[i] = p.label(valid.df$lmm.pval[i], valid.df$glmm.pval[i], nrow(results.df))
}
valid.df$p.cat = factor(valid.df$p.cat, levels=c("**/**", "**/*", "**/-", "*/*", "*/-", "-/-"))
g.mt.glmm = ggplot(valid.df, aes(x=glmm.coef, y=log(-log(glmm.pval)), label=label, color=p.cat)) + 
  geom_point() + geom_text_repel(max_overlaps=50, size=3, lineheight=0.75) +
  theme_light() + labs(title="MT gene count PGLMM", x ="PGLMM coefficient", y = "log(-log(p))", color="p profile")
g.mt.lmm = ggplot(valid.df, aes(x=lmm.coef, y=log(-log(lmm.pval)), label=label, color=p.cat)) + 
  geom_point() + geom_text_repel(max.overlaps=50, size=3, lineheight=0.75) +
  theme_light() + labs(title="MT gene count PLMM", x ="PLMM coefficient", y = "log(-log(p))", color="p profile")
grid.arrange(g.mt.glmm, g.mt.lmm)

sf = 2
png("mt-test-mixed.png", width=600*sf, height=400*sf, res=72*sf)
grid.arrange(g.mt.glmm, g.mt.lmm)
dev.off()

#results.df[which(results.df$colname == "salt.tolerance"),]
#cor.plot[[68]]
#results.df[which(results.df$colname == "woodiness"),]
#cor.plot[[122]]
#results.df[which(results.df$colname == "plant.growth.form"),]
#cor.plot[[22]]
#results.df[which(results.df$colname == "parasiteOf"),]
#cor.plot[[8]]
#results.df[which(results.df$colname == "habitat"),]
#cor.plot[[165]]

cor.plot[[48]]

png("mt-test-mixed-specifics.png", width=600, height=400)
grid.arrange(cor.plot[[68]], cor.plot[[122]], cor.plot[[22]], cor.plot[[8]], cor.plot[[165]], nrow=2)
dev.off()
