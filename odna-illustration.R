library(ggplot2)
library(ggtree)
library(phangorn)
library(phytools)
library(ggtreeExtra)
library(ggrepel)

tree = read.newick("tree-for-traits-clean-pt.phy")
ggtree(tree, layout='circular', branch.length='none') 

# read Kostas' dataset
df = read.table("PTFull22.txt", sep="\t", header=T, stringsAsFactors = TRUE)

# the dataset is read as factors by default -- we'll awkwardly convert between data types
df$Scientific.Name = as.character(df$Scientific.Name)
df$Ancestry = as.character(df$Ancestry)

convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

tree$tip.label = convname(tree$tip.label)

ggtree(tree, layout='circular', branch.length='none', alpha=0.1) %<+% df +
  geom_tippoint(aes(colour=parasiteOf, size = countsPT))

pt.counts = data.frame(Label = NULL, Count = NULL)
for(label in tree$tip.label) {
  ref = which(df$Scientific.Name == label)
  if(length(ref) != 0) {
    pt.counts = rbind(pt.counts, data.frame(Label=label, Count = df$countsPT[ref], 
                                            Predictor=df$parasiteOf[ref],
                      TextLabel=ifelse(df$parasiteOf[ref]=="parasiteOf", df$Scientific.Name[ref], "")))
  } else {
    print("wtf")
  }
}

sf=4
png("pt-tree-fig.png", width=1000*sf, height=1000*sf, res=72*sf)
ggtree(tree, layout='circular', branch.length='none') %<+% df +
  geom_fruit(data = pt.counts, geom=geom_point, pwidth=0.1,
             mapping = aes(y=Label, x=Predictor, colour=Predictor), 
             stat="identity", size=2) +
  geom_fruit(data = pt.counts, geom=geom_bar, pwidth=1.3,
             mapping = aes(y=Label, x=Count, fill=Predictor), 
             orientation="y", stat="identity", size=4) +
  theme(legend.position="none")
dev.off()

###############

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
tree.tips = convname(tree$tip.label)
tree.labels = c(tree$tip.label, tree$node.label)

# read Kostas' dataset
df = read.table("MTFull22.txt", sep="\t", header=T, stringsAsFactors = TRUE)

# the dataset is read as factors by default -- we'll awkwardly convert between data types
df$Scientific.Name = as.character(df$Scientific.Name)
df$Ancestry = as.character(df$Ancestry)

convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

tree$tip.label = convname(tree$tip.label)

ggtree(tree, layout='circular', branch.length='none', alpha=0.1) %<+% df +
  geom_tippoint(aes(colour=parasiteOf, size = countsMT))

mt.counts = data.frame(Label = NULL, Count = NULL)
for(label in tree$tip.label) {
  ref = which(df$Scientific.Name == label)
  if(length(ref) != 0) {
    mt.counts = rbind(mt.counts, data.frame(Label=label, Count = df$countsMT[ref], 
                                            Predictor=df$parasiteOf[ref],
                                            TextLabel=ifelse(df$parasiteOf[ref]=="parasiteOf", df$Scientific.Name[ref], "")))
  } else {
    print("wtf")
  }
}

sf=4
png("mt-tree-fig-cover.png", width=1000*sf, height=1000*sf, res=72*sf)
ggtree(tree, layout='circular', branch.length='none', size=0.8) %<+% df +
  geom_fruit(data = mt.counts, geom=geom_point, pwidth=0.1,
             mapping = aes(y=Label, x=Predictor, colour=Predictor), 
             stat="identity", size=2) +
  geom_fruit(data = mt.counts, geom=geom_bar, pwidth=1.3,
             mapping = aes(y=Label, x=Count, fill=Predictor), 
             orientation="y", stat="identity", size=4) +
  theme(legend.position="none")
dev.off()

