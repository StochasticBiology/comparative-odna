source("mt-test.R")
source("mt-test-ncbi.R")
source("pt-test.R")
source("pt-test-ncbi.R")

sf = 2
png("odna-corrs-plm.png", width=800*sf, height=600*sf, res=72*sf)
grid.arrange(g.mt.gene.plm, g.mt.orf.plm, 
             g.pt.gene.plm, g.pt.orf.plm,
             nrow=2)
dev.off()

sf = 2
png("odna-corrs-pglm.png", width=800*sf, height=600*sf, res=72*sf)
grid.arrange(g.mt.gene.pglm, g.mt.orf.pglm, 
             g.pt.gene.pglm, g.pt.orf.pglm,
             nrow=2)
dev.off()

source("mt-test-no-metazoa.R")
source("mt-test-ncbi-no-metazoa.R")
sf = 2
png("odna-corrs-no-metazoa.png", width=800*sf, height=600*sf, res=72*sf)
grid.arrange(g.mt.gene.nomet.plm, g.mt.orf.nomet.plm, 
             g.mt.gene.nomet.pglm, g.mt.orf.nomet.pglm,
             nrow=2)
dev.off()
