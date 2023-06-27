library(ggplot2)
library(ggrepel)
library(gridExtra)

# read mtDNA results file
df.mt = read.csv("mac-results-df-mt-pglm.csv", stringsAsFactors = FALSE)

# produce annotation labels for plotting
df.mt$pglm.label = ""
df.mt$pglm.p.cat = ""
for(i in 1:nrow(df.mt)) {
  if(df.mt$pglm.pval[i] < 0.05) {
  df.mt$pglm.label[i] = paste(c(df.mt$colname[i], ":\n", df.mt$positive.label[i], "\n", df.mt$mrca.positive[i]), collapse="")
  }
  df.mt$pglm.p.cat[i] = ifelse(df.mt$pglm.pval[i] < 0.05/nrow(df.mt), "BF p < 0.05", ifelse(df.mt$pglm.pval[i] < 0.05, "p < 0.05", "NS"))
}
# construct plot
plot.mt.pglm = ggplot(df.mt, aes(x=pglm.coef, y=log(-log(pglm.pval)), label=pglm.label, color=pglm.p.cat)) + 
         geom_point() + geom_text_repel(size=2, max.overlaps=20) + theme_classic() + ylim(0,NA)

# produce annotation labels for plotting
df.mt$plm.label = ""
df.mt$plm.p.cat = ""
for(i in 1:nrow(df.mt)) {
  if(df.mt$plm.pval[i] < 0.05) {
    df.mt$plm.label[i] = paste(c(df.mt$colname[i], ":\n", df.mt$positive.label[i], "\n", df.mt$mrca.positive[i]), collapse="")
  }
  df.mt$plm.p.cat[i] = ifelse(df.mt$plm.pval[i] < 0.05/nrow(df.mt), "BF p < 0.05", ifelse(df.mt$plm.pval[i] < 0.05, "p < 0.05", "NS"))
}
# construct plot
plot.mt.plm = ggplot(df.mt, aes(x=plm.coef, y=log(-log(plm.pval)), label=plm.label, color=plm.p.cat)) + 
  geom_point() + geom_text_repel(size=2, max.overlaps=20) + theme_classic() + ylim(0,NA)

# output the most statistically significant findings
#sigs.mt = df.mt[df.mt$p.cat != "NS",]
#for(i in 1:nrow(sigs.mt)) {
#  print(paste(c(sigs.mt$colname[i], ": ", sigs.mt$positive.label[i],
#              " within ", sigs.mt$mrca.positive[i], " -- coef ", round(sigs.mt$coef[i], digits=3), 
#              ", p ", round(sigs.mt$pval[i], digits=3)), collapse=""))
#}

grid.arrange(plot.mt.plm, plot.mt.pglm, nrow=2)

# read mtDNA results file
df.pt = read.csv("mac-results-df-pt-pglm.csv", stringsAsFactors = FALSE)

# produce annotation labels for plotting
df.pt$pglm.label = ""
df.pt$pglm.p.cat = ""
for(i in 1:nrow(df.pt)) {
  if(df.pt$pglm.pval[i] < 0.05) {
    df.pt$pglm.label[i] = paste(c(df.pt$colname[i], ":\n", df.pt$positive.label[i], "\n", df.pt$mrca.positive[i]), collapse="")
  }
  df.pt$pglm.p.cat[i] = ifelse(df.pt$pglm.pval[i] < 0.05/nrow(df.pt), "BF p < 0.05", ifelse(df.pt$pglm.pval[i] < 0.05, "p < 0.05", "NS"))
}
# construct plot
plot.pt.pglm = ggplot(df.pt, aes(x=pglm.coef, y=log(-log(pglm.pval)), label=pglm.label, color=pglm.p.cat)) + 
  geom_point() + geom_text_repel(size=2, max.overlaps=20) + theme_classic() + ylim(0,NA)

# produce annotation labels for plotting
df.pt$plm.label = ""
df.pt$plm.p.cat = ""
for(i in 1:nrow(df.pt)) {
  if(df.pt$plm.pval[i] < 0.05) {
    df.pt$plm.label[i] = paste(c(df.pt$colname[i], ":\n", df.pt$positive.label[i], "\n", df.pt$mrca.positive[i]), collapse="")
  }
  df.pt$plm.p.cat[i] = ifelse(df.pt$plm.pval[i] < 0.05/nrow(df.pt), "BF p < 0.05", ifelse(df.pt$plm.pval[i] < 0.05, "p < 0.05", "NS"))
}
# construct plot
plot.pt.plm = ggplot(df.pt, aes(x=plm.coef, y=log(-log(plm.pval)), label=plm.label, color=plm.p.cat)) + 
  geom_point() + geom_text_repel(size=2, max.overlaps=20) + theme_classic() + ylim(0,NA)

# output the most statistically significant findings
#sigs.pt = df.pt[df.pt$p.cat != "NS",]
#for(i in 1:nrow(sigs.pt)) {
#  print(paste(c(sigs.pt$colname[i], ": ", sigs.pt$positive.label[i],
#              " within ", sigs.pt$mrca.positive[i], " -- coef ", round(sigs.pt$coef[i], digits=3), 
#              ", p ", round(sigs.pt$pval[i], digits=3)), collapse=""))
#}

grid.arrange(plot.pt.plm, plot.pt.pglm, nrow=2)

##########

# read mtDNA results file without metazoa
df.mtnm = read.csv("mac-results-df-mt-no-metazoa-pglm.csv", stringsAsFactors = FALSE)

# produce annotation labels for plotting
df.mtnm$label = ""
df.mtnm$p.cat = ""
for(i in 1:nrow(df.mtnm)) {
  if(df.mtnm$pval[i] < 0.05) {
    df.mtnm$label[i] = paste(c(df.mtnm$colname[i], ":\n", df.mtnm$positive.label[i], "\n", df.mtnm$mrca.positive[i]), collapse="")
  }
  df.mtnm$p.cat[i] = ifelse(df.mtnm$pval[i] < 0.05/nrow(df.mtnm), "BF p < 0.05", ifelse(df.mtnm$pval[i] < 0.05, "p < 0.05", "NS"))
}
# construct plot
plot.mtnm = ggplot(df.mtnm, aes(x=coef, y=log(-log(pval)), label=label, color=p.cat)) + 
  geom_point() + geom_text_repel(size=2, max.overlaps=20) + theme_classic() + ylim(0,NA)

# output the most statistically significant findings
sigs.mtnm = df.mtnm[df.mtnm$p.cat != "NS",]
for(i in 1:nrow(sigs.mtnm)) {
  print(paste(c(sigs.mtnm$colname[i], ": ", sigs.mtnm$positive.label[i],
                " within ", sigs.mtnm$mrca.positive[i], " -- coef ", round(sigs.mtnm$coef[i], digits=3), 
                ", p ", round(sigs.mtnm$pval[i], digits=3)), collapse=""))
}

sf=2
png("plot-odna-corr-no-metazoa-pglm.png", width=400*sf, height=400*sf, res=72*sf)
plot.mtnm
dev.off()
