########################################## using ntp ####################################################
library(CMScaller)
#To prepare the required files for the NTP algorithm.
#RNA-seq platform data with TPM or FPKM should be log2(expr+1) and scale standardization.
#RNA-seq platform data with voom should be scale standardization.
#AFFY platform data with RMA should be scale standardization.
ccle.exp  <- read.csv("ccle.exprCRClog2tpm.csv",header = T, row.names = 1,check.names = F)
templates <- read.table("templates_CCCRC.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T,quote = "")

# apply template to ccle.exp
ccle.exp  <-ccle.exp[unique(templates$probe),]
ccle.exp <- t(scale(t(ccle.exp), scale = T, center = T))
ntp.ccle <- ntp(emat      = ccle.exp,
                templates = templates,
                doPlot    = T,
                nPerm     = 1000,
                distance  = "cosine",
                nCores    = 8,
                seed      = 202095,
                verbose   = T)
write.table(ntp.ccle, "ntp prediction in CCLE.txt",sep = "\t",row.names = T,col.names = NA,quote = F)