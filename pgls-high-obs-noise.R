library(phytools)
library(phangorn)
library(nlme)
library(ggplot2)
library(gridExtra)

# explores the power of PGLS to find correlations under some awkward conditions
# (i) the variables x and y are not well-behaved: x is binary and y discrete, and the rate of decrease of y depends on x
# (ii) x is not observed perfectly: p01 is the probability of seeing a 1 when x=0 and p10 the probability of seeing a 0 when x=1
#      [note for things like parasitism I think p01=0]
# x starts at 0 and has some probability of irreversibly changing to 1 at each step on the tree
# model.correlated determines whether y does depend on x or not; 10 samples are run for each case
# other variables are tree size and average number of x transitions

# this code focuses on PGLS and looks at more extreme data occlusion: up to 95% of positive predictor values are now observed as negatives

# initialise results frame and experiment counter
pvals = obs.stats = data.frame()
expt = 0

# we'll only look at one tree size here -- this structure is from other versions of the code
for(tree.size in 2**c(8)) {
  # create balanced tree with 2^n nodes
  my.tree = stree(tree.size, "balanced")
  my.tree$node.label = as.character(1:my.tree$Nnode)
  tree.labels = c(my.tree$tip.label, my.tree$node.label)
  my.tree$edge.length = rep(1, nrow(my.tree$edge))
  plot(my.tree)

  # loop through our model parameters
  for(mean.events in c(2, 4, 8)) {              # average number of predictor 0->1 transitions
    for(model.correlated in c(0,5,10,20)) {     # strength of predictor-response effect
      for(p01 in c(0)) {                        # false positive observation rate (0 for us)
        for(p10 in c(0,0.5,0.9,0.95)) {         # false negative observation rate
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
            
            # plot if we want
            #g.1 = ggplot(df, aes(x=factor(obs.x), y=y)) + geom_jitter(width=0.25) + geom_boxplot(alpha = 0.2)
            #g.2 = ggplot(tip.df, aes(x=factor(obs.x), y=y)) + geom_jitter(width=0.25) + geom_boxplot(alpha = 0.2)
            #grid.arrange(g.1, g.2, nrow=1)

            ####### basic, naive correlation
            basic.test = cor.test(tip.df$obs.x, tip.df$y)
            basic.pval = basic.test$p.value

            ####### PGLS under a Brownian model for correlations
            if(all(tip.df$obs.x == 0) | all(tip.df$obs.x == 1) | all(tip.df$y == 0) |
               (length(which(tip.df$y[tip.df$obs.x==1] > 0)) < 2) |
               (length(which(tip.df$y[tip.df$obs.x==0] > 0)) < 2)) { 
              pgls.pval =  NA 
            } else {
              my.correlation = corBrownian(phy=my.tree, form=~label)
             
                pgls.mod = gls(y~obs.x, correlation=my.correlation, data=tip.df, method="ML", na.action=na.omit)
                gls.mod = gls(y~obs.x, data=tip.df, method="ML", na.action=na.omit)
                pgls.coef = coef(pgls.mod)[2]
                pgls.pval = anova(pgls.mod)$'p-value'[2]
          
            }

            # store p-values in a data frame
            pvals = rbind(pvals, data.frame(tree.size=tree.size, model.correlated=model.correlated, mean.events=mean.events, 
                                p01=p01, p10=p10, expt=expt, basic.pval=basic.pval, pgls.pval=pgls.pval))
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

# plot annotation dataframe
label.text = data.frame(label=c("FP","TN","TP","FN"), x=c(1,1, 2.1,2.1)-0.3, y=c(5,-5,5,-5))

# PGLS performance
g.1 = ggplot(pvals[pvals$p01==0,], aes(x=factor(model.correlated), y=log10(-log10(pgls.pval)), color=noise.label)) + 
  geom_boxplot()+ geom_hline(yintercept=log10(-log10(0.05)), color="#888888")+ 
  geom_vline(xintercept=1.5, color="#888888")+
  geom_text(data=label.text,aes(x=x,y=y,label=label),color="#888888",size=3) +
  facet_grid(mean.events~tree.size) +
  theme_light() + xlab("True effect") + ylab("log(-log(p)) PGLS") + labs(color="Obs error")

# naive correlation performance
g.2 = ggplot(pvals[pvals$p01==0,], aes(x=factor(model.correlated), y=log10(-log10(basic.pval)), color=noise.label)) + 
  geom_boxplot()+ geom_hline(yintercept=log10(-log10(0.05)), color="#888888")+ 
  geom_vline(xintercept=1.5, color="#888888")+
  geom_text(data=label.text,aes(x=x,y=y,label=label),color="#888888",size=3) +
  facet_grid(mean.events~tree.size) + 
  theme_light() + xlab("True effect") + ylab("log(-log(p)) naive") + labs(color="Obs error")

# plot both
grid.arrange(g.1, g.2, nrow=1)

sf=2
png("pgls-high-obs-noise.png", width=600*sf, height=480*sf, res=72*sf)
grid.arrange(g.1, g.2, nrow=1)
dev.off()

# plot a single subplot for the main text
label.text = data.frame(label=c("FP","TN","TP","FN"), x=c(1,1, 2.1,2.1)-0.3, y=c(2.5,-2,2.5,-2))
g.sub = ggplot(pvals[pvals$p01==0 & pvals$mean.events==8,], aes(x=factor(model.correlated), y=log10(-log10(pgls.pval)), color=noise.label)) + 
  geom_boxplot()+ geom_hline(yintercept=log10(-log10(0.05)), color="#888888")+ 
  geom_vline(xintercept=1.5, color="#888888")+
  geom_text(data=label.text,aes(x=x,y=y,label=label),color="#888888",size=3) + 
  theme_light() + xlab("True effect") + ylab("log(-log(p)) naive") + labs(color="Obs error")

sf=2
png("pgls-high-obs-noise-subpanel.png", width=300*sf, height=300*sf, res=72*sf)
g.sub
dev.off()

# summarise info on observation statistics
obs.stats$mean.events = pvals$mean.events
obs.stats$tree.size = pvals$tree.size

ggplot(obs.stats, aes(x = mean.events, y=obs1/(obs0+obs1))) + geom_violin() + facet_grid(mean.events ~ tree.size)

png("pgls-high-obs-noise-stats.png", width=600, height=600)
ggplot(obs.stats, aes(x = mean.events, y=obs1/(obs0+obs1))) + geom_violin() + facet_grid(mean.events ~ tree.size)
dev.off()

