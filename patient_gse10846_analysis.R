# Author: Alex Bott
# Function: Analyze MPC1 and MPC2 expression in subtypes of DLBCL
# Date: 05/08/21

# https://www.biostars.org/p/221440/

# imports
library(GEOquery)
library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(hgu133plus2.db)

# classifications from Caro 2012
classifications <- read.csv('data/patient/id_subtype.csv')

# get data from GEO
gse <- getGEO('GSE10846')
gse <- gse[[1]]
se_gse <- makeSummarizedExperimentFromExpressionSet(gse)

# ensure everything is on the same platform
table(colData(se_gse)$platform_id == 'GPL570') 

# make a dataframe for all patients, then merge with classifications from Caro 2012
gsm <- rownames(colData(se_gse))
subtype <- as.data.frame(gsm)
colnames(subtype) <- 'id'
subtype <- left_join(subtype, classifications, by = 'id')

# add subtype to colData 
se_gse$subtype <- subtype
se_gse$subtype <- as.factor(subtype[ , 2])
desired_subtype <- c("BCR", "OxPhos")
cleaned <- se_gse[ , se_gse$subtype %in% desired_subtype]
cleaned$subtype <- droplevels(cleaned$subtype)

# UMAP plot following GEO2R
library(umap)
treatment <- se_gse$characteristics_ch1.9
ex <- assay(se_gse)
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, col = as.factor(treatment),  
     main="UMAP plot", xlab="", ylab="", pch=20, cex=1.5)

treatment2 <- cleaned$characteristics_ch1.9
ex2 <- assay(cleaned)
ump2 <- umap(t(ex2), n_neighbors = 15, random_state = 123)
plot(ump2$layout, col = as.factor(treatment2),  
     main="UMAP plot", xlab="", ylab="", pch=20, cex=1.5)

# the separating factor between the groups seemed to have been the drug regimen
# all the classified samples we are interested in are part of the R-CHOP-Like Regimen

# limma for differential expression
design <- model.matrix(~ 0 + cleaned$subtype)
colnames(design) <- c("BCR", "OxPhos")

fit <- lmFit(assay(cleaned), design)
cont.matrix <- makeContrasts(OxPhosVSBCR = OxPhos - BCR, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

output <- topTable(fit2, adjust="BH", n = nrow(fit2)) %>% 
  mutate(probe = rownames(.))

probes <- row.names(output)
decoder <- select(hgu133plus2.db, 
                  keys = probes, 
                  columns = c("SYMBOL", "GENENAME"),
                  keytype = "PROBEID")

output <- output %>% 
  left_join(decoder, by = c("probe" = "PROBEID"))

write.csv(output, "results/output/gse10846_limma_output.csv")

# probes for MPC1 and MPC2
mpc1_probe <- rowData(se_gse)$ID[which(rowData(se_gse)$Gene.Symbol == "MPC1")]
mpc2_probe <- rowData(se_gse)$ID[which(rowData(se_gse)$Gene.Symbol == "MPC2")]

# check mpc probes at https://gemma.msl.ubc.ca/arrays/showArrayDesign.html?id=4
# 218024_at aligns with CDS, 240362_at appears to be intronic
# 202427_s_at has lower percent identity when blasted, 227669_at aligns to 3' end and
# has 100% identity when blasted

# MPC1 
mpc1 <- subtype
mpc1$probe <- t(assay(se_gse['218024_at']))

# MCP2
mpc2 <- subtype
mpc2$probe <- t(assay(se_gse['227669_at']))

# plotting
ggplot(data = subset(mpc1, !is.na(subtype) & 
                       !subtype == 'Unassigned' & !subtype == 'HR'), 
       aes(x = subtype, y = probe)) + 
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(width = 0.25, alpha = 0.35) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  labs(title = 'MPC1', x = '', y = expression('log'[2]* ' signal intensity')) +
  theme_classic()

ggplot(data = subset(mpc2, !is.na(subtype) & 
                       !subtype == 'Unassigned' & !subtype == 'HR'), 
       aes(x = subtype, y = probe)) + 
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(width = 0.25, alpha = 0.35) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  labs(title = 'MPC2', x = '', y = 'Relative expression ' ~(log[2])) +
  theme_classic()

# extract median values by subtype
subset(mpc2, !is.na(subtype) & 
         !subtype == 'Unassigned' & !subtype == 'HR') %>% 
  group_by(subtype) %>% 
  summarize(med = median(probe))

# caro 2012 comparison
# PDHB - 211023_at, ACAT1 - 205412_at
subtype %>% 
  mutate(probe = t(assay(se_gse['211023_at']))) %>% 
  filter(!is.na(subtype) & !subtype == 'Unassigned' & !subtype == 'HR') %>% 
  ggplot(aes(x = subtype, y = probe)) + 
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(width = 0.25, alpha = 0.35) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  labs(title = 'PDHB', x = '', y = 'Relative expression ' ~(log[2])) +
  theme_classic()

subtype %>% 
  mutate(probe = t(assay(se_gse['205412_at']))) %>% 
  filter(!is.na(subtype) & !subtype == 'Unassigned' & !subtype == 'HR') %>% 
  ggplot(aes(x = subtype, y = probe)) + 
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(width = 0.25, alpha = 0.35) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  labs(title = 'ACAT1', x = '', y = 'Relative expression ' ~(log[2])) +
  theme_classic()

# plot
subtype %>% 
  mutate(probe = t(assay(se_gse['210337_s_at']))) %>% 
  filter(!is.na(subtype) & !subtype == 'Unassigned' & !subtype == 'HR') %>% 
  ggplot(aes(x = subtype, y = probe)) + 
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(width = 0.25, alpha = 0.35) +
  theme(axis.title.x = element_text(size = 13)) +
  theme(axis.title.y = element_text(size = 13)) +
  theme(axis.text.x = element_text(color = 'black')) +
  theme(axis.text.y = element_text(color = 'black')) +
  labs(title = 'ACLY', x = '', y = 'Relative expression ' ~(log[2])) +
  theme_classic()

# write MPC1 and MPC2 values
write.csv(subset(mpc1, !is.na(subtype) & 
         !subtype == 'Unassigned' & !subtype == 'HR'),
         "results/output/mpc1_microarray.csv",
         row.names = FALSE)
write.csv(subset(mpc2, !is.na(subtype) & 
         !subtype == 'Unassigned' & !subtype == 'HR'),
         "results/output/mpc2_microarray.csv",
         row.names = FALSE)

