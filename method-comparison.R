library(phytools)
library(phangorn)
library(nlme)
library(ggplot2)
library(gridExtra)

# explores the power of various comparative approaches to find correlations under some awkward conditions
# (i) the variables x and y are not well-behaved: x is binary and y discrete, and the rate of decrease of y depends on x
# (ii) x is not observed perfectly: p01 is the probability of seeing a 1 when x=0 and p10 the probability of seeing a 0 when x=1
#      [note for things like parasitism I think p01=0]
# x starts at 0 and has some probability of irreversibly changing to 1 at each step on the tree
# model.correlated determines whether y does depend on x or not; 10 samples are run for each case
# other variables are tree size and average number of x transitions

# initialise results frame and experiment counter
pvals = obs.stats = data.frame()
expt = 0

# loop over tree sizes
for(tree.size in 2**c(4, 5, 6)) {
  # create balanced tree with 2^n nodes
  my.tree = stree(tree.size, "balanced")
  my.tree$node.label = as.character(1:my.tree$Nnode)
  tree.labels = c(my.tree$tip.label, my.tree$node.label)
  my.tree$edge.length = rep(1, nrow(my.tree$edge))
  plot(my.tree)

  for(mean.events in c(2, 4, 8)) {              # average number of predictor 0->1 transitions
    for(model.correlated in c(0,5,10,20)) {     # strength of predictor-response effect
      for(p01 in c(0)) {                        # false positive observation rate (0 for us)
        for(p10 in c(0,0.1,0.25,0.5)) {         # false negative observation rate
          cat(paste(c(tree.size, mean.events, model.correlated, p01, p10, "\n"), collapse=" "))
          for(this.expt in 1:20) {
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
            
            # plot if desired
            #g.1 = ggplot(df, aes(x=factor(obs.x), y=y)) + geom_jitter(width=0.25) + geom_boxplot(alpha = 0.2)
            #g.2 = ggplot(tip.df, aes(x=factor(obs.x), y=y)) + geom_jitter(width=0.25) + geom_boxplot(alpha = 0.2)
            #grid.arrange(g.1, g.2, nrow=1)

            ####### basic, naive correlation
            basic.test = cor.test(tip.df$obs.x, tip.df$y)
            basic.pval = basic.test$p.value

            if(all(tip.df$obs.x == 0) | all(tip.df$obs.x == 1) | all(tip.df$y == 0) |
               (length(which(tip.df$y[tip.df$obs.x==1] > 0)) < 2) |
               (length(which(tip.df$y[tip.df$obs.x==0] > 0)) < 2)) { 
              pgls.pval =  pglm.pval  = NA
            } else {
            ####### PGLS under a Brownian model for correlations
            my.correlation = corBrownian(phy=my.tree, form=~label)
          
              pgls.mod = gls(y~obs.x, correlation=my.correlation, data=tip.df, method="ML", na.action=na.omit)
              gls.mod = gls(y~obs.x, data=tip.df, method="ML", na.action=na.omit)
              pgls.coef = coef(pgls.mod)[2]
              pgls.pval = anova(pgls.mod)$'p-value'[2]
              pgls.pval
          
            ####### PGLM with branch lengths
              rownames(tip.df) = tip.df$label
              pglm = phyloglm(y~obs.x, tip.df, my.tree, method="poisson_GEE")
              pglm.coef = summary(pglm)$coefficients[2,1]
              pglm.pval = summary(pglm)$coefficients[2,4]
              
              plm = phylolm(y~obs.x, tip.df, my.tree, model="BM")
              plm.coef = summary(plm)$coefficients[2,1]
              plm.pval = summary(plm)$coefficients[2,4]
            }
            
            ####### a non-parametric approach summarising response differences between close relatives (sisters, cousins, etc) with different predictor values

            #### this code labels each node in the tree with a trait value and a gene count
            # this proceeds as follows:
            # 1. label all nodes, and node descendants, from our positives list with positive values
            # 2. go back up the tree labelling positive any parents with all positive children
            # 3. label all tips with their gene counts
            # 4. go back up the tree, labelling parent nodes with the average of their childrens' gene counts

            # 1.
            # set all nodes in the tree to negative value of trait by default
            tree.vals = rep(0, length(tree.labels))

            #cat("------- Painting nodes and descendants...\n  ")

            positives.list = tip.df$label[tip.df$obs.x == 1]
            if(length(positives.list) == 0) {
              wilcox.pval = NA
              boot.pval = NA
            } else {
              # go through list of positives, painting positive trait values onto the tree
              for(i in 1:length(positives.list)) {
                name = tolower(positives.list[i])
                if(name != "") {
                  # find nodes with this positive label
                  tree.refs = which(tree.labels == name)
                  for(tree.ref in tree.refs) {
                    # set this node's trait value to positive
                    tree.vals[tree.ref] = 1
                    #cat(tree.labels[tree.ref])
                    #cat("\n")
                    # set descendant nodes' trait values to positive
                    descendant.refs = getDescendants(my.tree, tree.ref)
                    #cat(paste(c("    ", tree.labels[descendant.refs]), collapse=", "))
                    #cat("\n  ")
                    tree.vals[descendant.refs] = 1
                  }
                }
              }

              #cat("\n------- Painting ancestors...\n  ")

              # 2.
              # this goes back up the tree, painting a positive value on any parent node will all positive descendants
              # ideally this shouldn't be necessary, but there are instances where all species are labelled positive without their clade being positive
              # this snippet is very inefficient and can be sped up several ways
              change = T
              # while we're still not converged
              while(change == T) {
                change = F
                # loop through all nodes
                for(tree.ref in 1:length(tree.labels)) {
                  if(tree.vals[tree.ref] == 0) {
                    # if this node is negative, check to see if its descendants are all positive
                    descendant.refs = getDescendants(my.tree, tree.ref)
                    if(all(tree.vals[descendant.refs] == 1)) {
                      # if so, make the parent positive
                      #cat(tree.labels[tree.ref])
                      #cat("\n  ")
                      tree.vals[tree.ref] = 1
                      change = T
                    }
                  }
                }
              }

              # 3.
              # now get response variable for tips
              tree.y = rep(-999, length(tree.labels))
              # loop through tips
              for(tip.ref in 1:length(my.tree$tip.label)) {
                # grab barcode entry corresponding to this node's label
                tree.y[tip.ref] = df$y[df$label == my.tree$tip.label[tip.ref]]
              }

              # 4.
              # now move back up the tree, labelling parent nodes with the average of their childrens' gene count
              change = T
              # while we're not converged
              while(change == T) {
                change = F
                # go through nodes
                for(tree.ref in 1:length(tree.labels)) {
                  if(tree.y[tree.ref] == -999) {
                    # this node hasn't got a count label yet, so if its children all have one, use their mean
                    descendant.refs = Children(my.tree, tree.ref)
                    if(all(tree.y[descendant.refs] != -999)) {
                      tree.y[tree.ref] = mean(tree.y[descendant.refs])
                      change = T
                    }
                  }
                }
              }


              #### now we can find pairs of sibling nodes whose trait value differs

              #cat("\n------- Identifying sibling pairs...\n  ")

              sort.df = df[order(df$ref),]

              # initialise a dataframe and "done" marker 
              pairs = data.frame(ego=NULL,myval=NULL,theirval=NULL,mycount=NULL,theircount=NULL)
              included = rep(0, length(tree.labels))

              # loop through all nodes
              for(tree.ref in 1:length(tree.labels)) {
                # get this node's set of siblings (including self)
                siblings = Siblings(my.tree, tree.ref, include.self=T)
  
                # sets of siblings with opposite trait value and with same trait value
                opposites = which(sort.df$obs.x[siblings] != sort.df$obs.x[tree.ref])
                sames = which(sort.df$obs.x[siblings] == sort.df$obs.x[tree.ref])
  
                # if we have elements in both sets, which include no nodes that have yet been considered
                if(length(opposites) > 0 & all(included[siblings[sames]] == 0) & all(included[siblings[opposites]] == 0)) {
                  # add the pair of the sets' mean gene counts to our list of comparisons
                  #cat(paste(c("{", paste(c(tree.labels[siblings[sames]]), collapse=", "), "} vs {", paste(c(tree.labels[siblings[opposites]]), collapse=", ")  , 
                  #            "}\n     ==> {",  
                  #            c(sort.df$y[siblings[sames]]), " } vs {", c(sort.df$y[siblings[opposites]]), "}"), collapse=" "))
                  #cat("\n  ")
                  pairs = rbind(pairs, data.frame(ego=tree.ref, 
                                                  myval=mean(sort.df$obs.x[siblings[sames]]), theirval=mean(sort.df$obs.x[siblings[opposites]]), 
                                                  mycount=mean(sort.df$y[siblings[sames]]), theircount=mean(sort.df$y[siblings[opposites]])))
                  # note that we've included all these siblings.
                  # is this fair?
                  included[siblings] = 1
                }
              }

              #cat("\n------- Tests...\n  ")

              # organise dataframe so that we're looking at (without trait)-(with trait) differences
              pairs$diff = NULL
              if(nrow(pairs) == 0) { 
                wilcox.pval = boot.pval = NA
              } else {
                for(i in 1:nrow(pairs)) {
                  if(pairs$myval[i] == 1) {
                    pairs$diff[i] = pairs$mycount[i]-pairs$theircount[i]
                  } else {
                    pairs$diff[i] = pairs$theircount[i]-pairs$mycount[i]
                  }
                }
                
                #cat(paste(c(length(pairs$diff), " pair comparisons\n(Positive-negative):\n  Mean ", mean(pairs$diff), "\n  SD ", sd(pairs$diff)), collapse=""))

                wilcox = wilcox.test(pairs$diff)
                wilcox.pval = wilcox$p.value

                # bootstrap for mean difference
                means = NULL
                nboot = 10000
                for(boot in 1:nboot) {
                  booted = sample(pairs$diff, length(pairs$diff), replace=T)
                  means = c(means, mean(booted))
                }
                if(mean(pairs$diff) > 0) { samples = sum(means < 0)/nboot } else { samples = sum(means > 0)/nboot }
                #cat(paste(c("Boot samples past zero = ", samples, " of ", nboot, "\n"), collapse=""))
                boot.pval = samples
                if(boot.pval == 0) { boot.pval = 1/nboot }
              }
            }

            # store p-values in a data frame
            pvals = rbind(pvals, data.frame(tree.size=tree.size, model.correlated=model.correlated, mean.events=mean.events, 
                                p01=p01, p10=p10, expt=expt, basic.pval=basic.pval, pgls.pval=pgls.pval, pglm.pval=pglm.pval, 
                                plm.pval=plm.pval, wilcox.pval=wilcox.pval, boot.pval=boot.pval))
          }
        }
      }
    }
  }
}

# for plotting, produce a label describing the observation noise for each reading
pvals$noise.label = ""
for(i in 1:nrow(pvals)) {
  #pvals$noise.label[i] = paste(c(pvals$p01[i], pvals$p10[i]), collapse=" ")
  pvals$noise.label[i] = paste(c(pvals$p10[i]), collapse=" ")
}

# plot annotation data frame
label.text = data.frame(label=c("FP","TN","TP","FN"), x=c(1,1, 2.1,2.1)-0.4, y=c(1,-1,1,-1))

# PGLS results
g.1 = ggplot(pvals[pvals$p01==0,], aes(x=factor(model.correlated), y=log10(-log10(pgls.pval)), color=noise.label)) + 
  geom_boxplot() + geom_hline(yintercept=log10(-log10(0.05)), color="#888888")+ 
  geom_vline(xintercept=1.5, color="#888888")+
  geom_text(data=label.text,aes(x=x,y=y,label=label),color="#888888",size=3) +
  facet_grid(mean.events~tree.size) + 
  theme_light() + xlab("True effect") + ylab("log(-log(p)) PGLS") + labs(color="Obs error")

# naive correlation results
g.2 = ggplot(pvals[pvals$p01==0,], aes(x=factor(model.correlated), y=log10(-log10(basic.pval)), color=noise.label)) + 
  geom_boxplot() + geom_hline(yintercept=log10(-log10(0.05)), color="#888888")+ 
  geom_vline(xintercept=1.5, color="#888888")+
  geom_text(data=label.text,aes(x=x,y=y,label=label),color="#888888",size=3) +
  facet_grid(mean.events~tree.size) +
  theme_light() + xlab("True effect") + ylab("log(-log(p)) naive") + labs(color="Obs error")

# relative comparison via Wilcoxon
g.3 = ggplot(pvals[pvals$p01==0,], aes(x=factor(model.correlated), y=log10(-log10(wilcox.pval)), color=noise.label)) + 
  geom_boxplot() + geom_hline(yintercept=log10(-log10(0.05)), color="#888888")+ 
  geom_vline(xintercept=1.5, color="#888888")+
  geom_text(data=label.text,aes(x=x,y=y,label=label),color="#888888",size=3) +
  facet_grid(mean.events~tree.size) + 
  theme_light() + xlab("True effect") + ylab("log(-log(p)) relatives-Wilcoxon") + labs(color="Obs error")

# relative comparison via bootstrap
g.4 = ggplot(pvals[pvals$p01==0,], aes(x=factor(model.correlated), y=log10(-log10(boot.pval)), color=noise.label)) + 
  geom_boxplot() + geom_hline(yintercept=log10(-log10(0.05)), color="#888888")+ 
  geom_vline(xintercept=1.5, color="#888888")+
  geom_text(data=label.text,aes(x=x,y=y,label=label),color="#888888",size=3) +
  facet_grid(mean.events~tree.size) + 
  theme_light() + xlab("True effect") + ylab("log(-log(p)) relatives-bootstrap") + labs(color="Obs error")

# PGLM
g.5 = ggplot(pvals[pvals$p01==0,], aes(x=factor(model.correlated), y=log10(-log10(pglm.pval)), color=noise.label)) + 
  geom_boxplot() + geom_hline(yintercept=log10(-log10(0.05)), color="#888888")+ 
  geom_vline(xintercept=1.5, color="#888888")+
  geom_text(data=label.text,aes(x=x,y=y,label=label),color="#888888",size=3) +
  facet_grid(mean.events~tree.size) + 
  theme_light() + xlab("True effect") + ylab("log(-log(p)) PGLM") + labs(color="Obs error")

# PLM
g.6 = ggplot(pvals[pvals$p01==0,], aes(x=factor(model.correlated), y=log10(-log10(plm.pval)), color=noise.label)) + 
  geom_boxplot() + geom_hline(yintercept=log10(-log10(0.05)), color="#888888")+ 
  geom_vline(xintercept=1.5, color="#888888")+
  geom_text(data=label.text,aes(x=x,y=y,label=label),color="#888888",size=3) +
  facet_grid(mean.events~tree.size) + 
  theme_light() + xlab("True effect") + ylab("log(-log(p)) PLM") + labs(color="Obs error")

# plot all together
grid.arrange(g.1, g.6, g.3, g.4, g.5, g.2, nrow=3)

sf = 2
png("method-comparison.png", width=1000*sf, height=1000*sf, res=72*sf)
grid.arrange(g.1, g.6, g.3, g.4, g.5, g.2, nrow=3)
dev.off()

# summarise info on observation statistics
obs.stats$mean.events = pvals$mean.events
obs.stats$tree.size = pvals$tree.size

ggplot(obs.stats, aes(x = mean.events, y=obs1/(obs0+obs1))) + geom_violin() + facet_grid(mean.events ~ tree.size)

png("method-comparison-obs.png", width=600, height=600)
ggplot(obs.stats, aes(x = mean.events, y=obs1/(obs0+obs1))) + geom_violin() + facet_grid(mean.events ~ tree.size)
dev.off()

