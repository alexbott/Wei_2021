# Author: Alex Bott
# Function: Analyze MPC1 and MPC2 expression in subtypes of DLBCL
# Date: 05/08/21

# imports
library(GEOquery)
library(SummarizedExperiment)
library(tidyverse)

# classifications from Caro 2012
classifications <- read.csv('data/patient/id_subtype.csv')

# get data from GEO
gse <- getGEO('GSE10846')
gse <- gse[[1]]
se_gse <- makeSummarizedExperimentFromExpressionSet(gse)

# ensure everything is on the same platform
table(colData(se_gse)$platform_id == 'GPL570') 

# probes for MPC1 and MPC2
mpc1_probe <- rowData(se_gse)$ID[which(rowData(se_gse)$Gene.Symbol == "MPC1")]
mpc2_probe <- rowData(se_gse)$ID[which(rowData(se_gse)$Gene.Symbol == "MPC2")]

# check coefficient of variance on various probes
cov <- function(x) {
  sd(x) / mean(x)
  }

apply(assay(se_gse[mpc1_probe]), 1, cov)
apply(assay(se_gse[mpc2_probe]), 1, cov)

# make a dataframe for all patients, then merge with classifications from Caro 2012
gsm <- rownames(colData(se_gse))
subtype <- as.data.frame(gsm)
colnames(subtype) <- 'id'
subtype <- left_join(subtype, classifications, by = 'id')

# add subtype to colData 
se_gse$subtype <- subtype
se_gse$subtype <- as.factor(subtype[ , 2])

# MPC1 
mpc1 <- subtype
mpc1$p1 <- t(assay(se_gse['218024_at']))
mpc1$p2 <- t(assay(se_gse['240362_at']))

# MCP2
mpc2 <- subtype
mpc2$p1 <- t(assay(se_gse['202427_s_at']))
mpc2$p2 <- t(assay(se_gse['227669_at']))


# plotting
ggplot(data = subset(mpc1, !is.na(subtype) & 
                       !subtype == 'Unassigned' & !subtype == 'HR'), 
       aes(x = subtype, y = p1)) + 
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(width = 0.25, alpha = 0.35) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  labs(title = 'MPC1', x = '', y = 'Relative expression ' ~(log[2])) +
  theme_classic()
