library(phytools)
library(phangorn)
library(nlme)
library(ggplot2)
library(gridExtra)

library(ggtree)

my.seed = 4

set.seed(my.seed)

# create and plot example trees and data to illustrate the synthetic evolutionary work

# initialise results frame and experiment counter
pvals = obs.stats = data.frame()
expt = 0

tree.size = 32
mean.events = 4
model.correlated = 0
p01 = 0
p10 = 0.5

  # create balanced tree with 2^n nodes
  my.tree = stree(tree.size, "balanced")
  my.tree$node.label = as.character(1:my.tree$Nnode)
  tree.labels = c(my.tree$tip.label, my.tree$node.label)
  my.tree$edge.length = rep(1, nrow(my.tree$edge))
  plot(my.tree)

 
            # get root reference, and add it to a "to-do" list for simulating traits
            my.root = getRoot(my.tree)
            to.do = c(my.root)
            # initialise state list
            df = data.frame(ref=my.root, label=tree.labels[my.root], x=0, y=100)
            # while we still have vertices to simulate
            while(length(to.do) > 0) {
              # initialise a new to-do list for next iteration
              new.to.do = c()
              # loop through current to-do list
              for(i in to.do) {
                # get state of this vertex
                x = df$x[df$ref == i]
                y = df$y[df$ref == i]
                # loop over its children
                for(j in Children(my.tree, i)) {
                  # switch x state with some probability
                  if(runif(1) < mean.events/tree.size | x == 1) { new.x = 1 } else { new.x = 0 }
                  # reduce y state by an x-dependent amount
                  if(new.x == 1) { dy = round(runif(1, min=-5-model.correlated, max=-model.correlated)) } 
                  else { dy = round(runif(1, min=-5, max=0)) } 
                  new.y = y + dy
                  
                  # add this child to to state list, and to next iteration's to-do
                  df = rbind(df, data.frame(ref=j, label=tree.labels[j], x=new.x, y=new.y))
                  new.to.do = c(new.to.do, j)
                }
              }
              # update to-do list
              to.do = new.to.do
            }

            # let's use obs.x to model our understanding of the predictor variable
            # this imposes false-positive and false-negative observation probabilities
            df$obs.x = df$x
            obs.0 = which(df$obs.x == 0)
            obs.1 = which(df$obs.x == 1)
            obs.01 = sample(obs.0, p01*length(obs.0))
            obs.10 = sample(obs.1, p10*length(obs.1))
            df$obs.x[obs.01] = 1
            df$obs.x[obs.10] = 0

            # pull just the tip observations from the tree
            tip.df = df[df$ref <= length(my.tree$tip.label),]
            tip.df$y[tip.df$y < 0] = 0
            # record observation stats (how many real/observed positives/negatives)
            obs.stats = rbind(obs.stats, data.frame(expt=expt, true0=length(which(tip.df$x == 0)),
                                                    true1=length(which(tip.df$x == 1)),
                                                    obs0=length(which(tip.df$obs.x == 0)),
                                                    obs1=length(which(tip.df$obs.x == 1))))
  

plot.tip.df = tip.df
plot.tip.df$ref = NULL
colnames(plot.tip.df) = c("df.label", "df.x", "df.y", "df.obs.x")
g.tree.1 = ggtree(my.tree) %<+% plot.tip.df + 
   geom_tippoint(aes(shape = factor(df.x), color = factor(df.x), size = df.y), position=position_nudge(x=0.25)) +
   geom_tippoint(aes(shape = factor(df.obs.x), color=factor(df.obs.x), size = df.y), position=position_nudge(x=0.5)) +
  theme(legend.position="none")

 ###########
 
set.seed(my.seed)

tree.size = 32
death = 1
mean.events = 4
model.correlated = 10
p01 = 0
p10 = 0.5
               # create random phylogeny with 2^n nodes from birth-death process parameterised as above
               my.tree = rphylo(tree.size, birth=1, death=death)
               my.tree$node.label = as.character(1:my.tree$Nnode)
               tree.labels = c(my.tree$tip.label, my.tree$node.label)
               my.tree$edge.length = my.tree$edge.length/mean(my.tree$edge.length)
               
               expt = expt+1
               # get root reference, and add it to a "to-do" list for simulating traits
               my.root = getRoot(my.tree)
               to.do = c(my.root)
               # initialise state list
               df = data.frame(ref=my.root, label=tree.labels[my.root], x=0, y=100)
               # while we still have vertices to simulate
               while(length(to.do) > 0) {
                 # initialise a new to-do list for next iteration
                 new.to.do = c()
                 # loop through current to-do list
                 for(i in to.do) {
                   # get state of this vertex
                   x = df$x[df$ref == i]
                   y = df$y[df$ref == i]
                   this.outgoing.edges = which(my.tree$edge[,1] == i)
                   # loop over its children
                   for(j in this.outgoing.edges) {
                     this.branch.length = my.tree$edge.length[j]
                     this.child = my.tree$edge[j,2]
                     # switch x state with some probability
                     #print(c("Processing edge ", j, " child ", this.child))
                     if(runif(1) < mean.events/tree.size | x == 1) { new.x = 1 } else { new.x = 0 }
                     # reduce y state by an x- and time-dependent amount
                     if(new.x == 1) { dy = round(this.branch.length*runif(1, min=-5-model.correlated, max=-model.correlated)) } 
                     else { dy = round(this.branch.length*runif(1, min=-5, max=0)) } 
                     new.y = y + dy
                     
                     # add this child to to state list, and to next iteration's to-do
                     df = rbind(df, data.frame(ref=this.child, label=tree.labels[this.child], x=new.x, y=new.y))
                     new.to.do = c(new.to.do, this.child)
                   }
                 }
                 # update to-do list
                 to.do = new.to.do
               }
               
               # let's use obs.x to model our understanding of the predictor variable
               # this imposes false-positive and false-negative observation probabilities
               df$obs.x = df$x
               obs.0 = which(df$obs.x == 0)
               obs.1 = which(df$obs.x == 1)
               obs.01 = sample(obs.0, p01*length(obs.0))
               obs.10 = sample(obs.1, p10*length(obs.1))
               df$obs.x[obs.01] = 1
               df$obs.x[obs.10] = 0
               
               # pull just the tip observations from the tree
               tip.df = df[df$ref <= length(my.tree$tip.label),]
               tip.df$label = as.character(tip.df$label)
               tip.df$y[tip.df$y < 0] = 0
               
               # record observation stats (how many real/observed positives/negatives)
               obs.stats = rbind(obs.stats, data.frame(expt=expt, true0=length(which(tip.df$x == 0)),
                                                       true1=length(which(tip.df$x == 1)),
                                                       obs0=length(which(tip.df$obs.x == 0)),
                                                       obs1=length(which(tip.df$obs.x == 1))))
               
               plot.tip.df = tip.df
               plot.tip.df$ref = NULL
               colnames(plot.tip.df) = c("df.label", "df.x", "df.y", "df.obs.x")
  g.tree.2 =              ggtree(my.tree) %<+% plot.tip.df + 
                 geom_tippoint(aes(shape = factor(df.x), color = factor(df.x), size = df.y), position=position_nudge(x=0.25)) +
                 geom_tippoint(aes(shape = factor(df.obs.x), color=factor(df.obs.x), size = df.y), position=position_nudge(x=0.5)) +
    theme(legend.position="none")
         
               ############
               set.seed(my.seed)
               tree.size = 32
               death = 0.1
               mean.events = 4
               model.correlated = 10
               p01 = 0
               p10 = 0.5
               # create random phylogeny with 2^n nodes from birth-death process parameterised as above
               my.tree = rphylo(tree.size, birth=1, death=death)
               my.tree$node.label = as.character(1:my.tree$Nnode)
               tree.labels = c(my.tree$tip.label, my.tree$node.label)
               my.tree$edge.length = my.tree$edge.length/mean(my.tree$edge.length)
               
               expt = expt+1
               # get root reference, and add it to a "to-do" list for simulating traits
               my.root = getRoot(my.tree)
               to.do = c(my.root)
               # initialise state list
               df = data.frame(ref=my.root, label=tree.labels[my.root], x=0, y=100)
               # while we still have vertices to simulate
               while(length(to.do) > 0) {
                 # initialise a new to-do list for next iteration
                 new.to.do = c()
                 # loop through current to-do list
                 for(i in to.do) {
                   # get state of this vertex
                   x = df$x[df$ref == i]
                   y = df$y[df$ref == i]
                   this.outgoing.edges = which(my.tree$edge[,1] == i)
                   # loop over its children
                   for(j in this.outgoing.edges) {
                     this.branch.length = my.tree$edge.length[j]
                     this.child = my.tree$edge[j,2]
                     # switch x state with some probability
                     #print(c("Processing edge ", j, " child ", this.child))
                     if(runif(1) < mean.events/tree.size | x == 1) { new.x = 1 } else { new.x = 0 }
                     # reduce y state by an x- and time-dependent amount
                     if(new.x == 1) { dy = round(this.branch.length*runif(1, min=-5-model.correlated, max=-model.correlated)) } 
                     else { dy = round(this.branch.length*runif(1, min=-5, max=0)) } 
                     new.y = y + dy
                     
                     # add this child to to state list, and to next iteration's to-do
                     df = rbind(df, data.frame(ref=this.child, label=tree.labels[this.child], x=new.x, y=new.y))
                     new.to.do = c(new.to.do, this.child)
                   }
                 }
                 # update to-do list
                 to.do = new.to.do
               }
               
               # let's use obs.x to model our understanding of the predictor variable
               # this imposes false-positive and false-negative observation probabilities
               df$obs.x = df$x
               obs.0 = which(df$obs.x == 0)
               obs.1 = which(df$obs.x == 1)
               obs.01 = sample(obs.0, p01*length(obs.0))
               obs.10 = sample(obs.1, p10*length(obs.1))
               df$obs.x[obs.01] = 1
               df$obs.x[obs.10] = 0
               
               # pull just the tip observations from the tree
               tip.df = df[df$ref <= length(my.tree$tip.label),]
               tip.df$label = as.character(tip.df$label)
               tip.df$y[tip.df$y < 0] = 0
               
               # record observation stats (how many real/observed positives/negatives)
               obs.stats = rbind(obs.stats, data.frame(expt=expt, true0=length(which(tip.df$x == 0)),
                                                       true1=length(which(tip.df$x == 1)),
                                                       obs0=length(which(tip.df$obs.x == 0)),
                                                       obs1=length(which(tip.df$obs.x == 1))))
               
               plot.tip.df = tip.df
               plot.tip.df$ref = NULL
               colnames(plot.tip.df) = c("df.label", "df.x", "df.y", "df.obs.x")
          g.tree.3=     ggtree(my.tree) %<+% plot.tip.df + 
                 geom_tippoint(aes(shape = factor(df.x), color = factor(df.x), size = df.y), position=position_nudge(x=0.25)) +
                 geom_tippoint(aes(shape = factor(df.obs.x), color=factor(df.obs.x), size = df.y), position=position_nudge(x=0.5)) +
            theme(legend.position="none")
          
          sf=1
          png("pgls-illustrations.png", width=1200*sf, height=500*sf, res=72*sf)
          grid.arrange(g.tree.1, g.tree.2, g.tree.3, nrow=1)
          dev.off()
               
               