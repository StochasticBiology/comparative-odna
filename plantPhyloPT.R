library(phytools)
library(phyreg)
library(ape)
library(geiger)
library(nlme)
library(caper)
library(stringr)
library(data.table)

library(taxizedb)
library(U.PhyloMaker)

convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}


df = read.table("PTFull22.txt", sep="\t", header=T)




megatree <- read.tree('https://raw.githubusercontent.com/megatrees/plant_20221117/main/plant_megatree.tre')
sp.list <- df[,1]

ids <- name2taxid(sp.list, out_type="summary")
sp.list=vector()
aa=classification(ids$id)

gen.list=data.frame(genus=character(),family=character())
for (i in 1:nrow(ids)) {
  if (length(aa[[i]][which(aa[[i]][2]=="kingdom"),1])>0) {
    if (aa[[i]][which(aa[[i]][2]=="kingdom"),1]=="Viridiplantae") { # class Mammalia
      gen.list[nrow(gen.list)+1,]=c(aa[[i]][which(aa[[i]][2]=="species"),1],aa[[i]][which(aa[[i]][2]=="family"),1])
    }
  }
}

sp.list=gen.list
sp.list[,1]=gen.list[,1]
sp.list[,2]=word(sp.list[,1], 1) 
sp.list[,3:5]=NA
colnames(sp.list) <- c("species","genus","family","species.relative","genus.relative")
gen.list[,1]=word(gen.list[,1], 1) 


result <- phylo.maker(sp.list=sp.list, tree=megatree, gen.list=gen.list, nodes.type = 1, scenario = 3)
result

# remove species not inserted 
bb=result$sp.list
length(which(bb$output.note!="no insertion"))
hits=bb[which(bb$output.note!="no insertion"),1]
truePhyloTree=result$phylo

df2=df
df2=df[match(tolower(hits),tolower(df$Scientific.Name)),]
df2=df2[which(!is.na(df2$Scientific.Name)),]


# clear some lables found in the tree but not in our specieslist
df=df2
tree=truePhyloTree
leafs=tree$tip.label
toDrop=setdiff(convname(tree$tip.label),df[,1])
toDrop=leafs[which(!is.na(match(convname(tree$tip.label),toDrop)))]
tree<-drop.tip(tree,toDrop)

write.table(df, "PTspeciesHits.txt", quote=FALSE, row.names=FALSE, sep = "\t",na="")
write.tree(tree,file = "pt-pruned-plantPhylo.phy")
