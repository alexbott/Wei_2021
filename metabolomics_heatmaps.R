# Author: Alex Bott
# Function: Analyze metabolomics data and generate heatmaps
# Date: 06/18/21

library(tidyverse)
library(pheatmap)
library(viridis)


# Normalized --------------------------------------------------------------

normRaw <- read.csv("data/metabolomics/20201005_PengWei_normalized.csv")
metabolites <- normRaw$Compound.Name
hmdb <- normRaw$HMDB
kegg <- normRaw$KEGG

# Heatmaps ----------------------------------------------------------------

# Each time point - heatmap for 4h, another for 8h, etc
# Vehice and UK - Suspension and matrigel
# Suspension UK has no growth phenotype, but matrigel and matrigel UK does
# Even though no growth phenotype, the impact of UK is similar between suspension and matrigel

heat <- normRaw
row.names(heat) <- heat$Compound.Name
heat <- heat[,-c(1:4)]

heat <- heat %>% 
  select(-contains("KO")) %>% 
  select(-contains("0hr"))

renameColsHeat <- paste(sub("^X[[:digit:]]+_", "", colnames(heat)),
                        sep = ".",
                        rep(seq(1:3), 12))
renameColsHeat <- sub("M", "Matrigel", renameColsHeat)
renameColsHeat <- sub("sus", "Suspension", renameColsHeat)
renameColsHeat <- sub("WT", "Veh", renameColsHeat)
colnames(heat) <- renameColsHeat


heat4 <- heat %>% 
  select(contains("4hr")) %>% 
  select(!contains("24hr"))
heat8 <- heat %>% 
  select(grep("8hr", colnames(heat)))
heat12 <- heat %>% 
  select(grep("12hr", colnames(heat)))
heat24 <- heat %>% 
  select(grep("24hr", colnames(heat)))

heatWTMat <- heat %>% 
  select(!contains("UK"))



annotation <- data.frame(Environment = c(rep("Matrigel", 6), rep("Suspension", 6)),
                         Treatment = c(rep("Vehicle", 3), rep("UK", 3), rep("Vehicle", 3), rep("UK", 3)))
annotation$Environment <- factor(annotation$Environment, levels = c("Suspension", "Matrigel"))
annotation$Treatment <- factor(annotation$Treatment, levels = c("Vehicle", "UK"))



row.names(annotation) <- colnames(heat24)
pheat24 <- heat24 %>% 
  select(-c("Veh.Matrigel.24hr.1")) %>% 
  log2() %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",
           show_rownames = FALSE,
           cluster_cols = F,
           cluster_rows = T,
           annotation = annotation,
           col = inferno(100))

rowOrder24 <- pheat24$tree_row$order

a <- heat24[rowOrder24, ]
b <- heat12[rowOrder24, ]
c <- heat8[rowOrder24, ]
d <- heat4[rowOrder24, ]


# save heatmaps
pdf(file = "results/figs/heatmap24h.pdf",
    width = 4, height = 4)

row.names(annotation) <- colnames(a)
aPlot <- a %>% 
  select(-c("Veh.Matrigel.24hr.1")) %>% 
  log2() %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",
           show_rownames = F,
           show_colnames = F,
           cluster_cols = T,
           cluster_rows = F,
           annotation = annotation,
           col = inferno(100))
dev.off()

pdf(file = "results/figs/heatmap12h.pdf",
    width = 4, height = 4)

row.names(annotation) <- colnames(b)
bPlot <- b %>% 
  log2() %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",
           show_rownames = F,
           show_colnames = F,
           cluster_cols = T,
           cluster_rows = F,
           annotation = annotation,
           col = inferno(100))
dev.off()

pdf(file = "results/figs/heatmap8h.pdf",
    width = 4, height = 4)

row.names(annotation) <- colnames(c)
cPlot <- c %>% 
  select(-c("Veh.Suspension.8hr.2")) %>% 
  log2() %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",
           show_rownames = F,
           show_colnames = F,
           cluster_cols = T,
           cluster_rows = F,
           annotation = annotation,
           col = inferno(100))
dev.off()

pdf(file = "results/figs/heatmap4h.pdf",
    width = 4, height = 4)

row.names(annotation) <- colnames(d)
dPlot <- d %>% 
  log2() %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",
           show_rownames = F,
           show_colnames = F,
           cluster_cols = T,
           cluster_rows = F,
           annotation = annotation,
           col = inferno(100))
dev.off()

# reverse 4h plot
library(seriation)
library(dendextend)

col_dend <- dPlot[[2]]
col_dend <- dendextend::rotate(col_dend, 
                               order = c(dPlot$tree_col$labels[rev(dPlot$tree_col$order)]))

pdf(file = "results/figs/heatmap4h_flipped.pdf",
    width = 4, height = 4)

row.names(annotation) <- colnames(d)
ePlot <- d %>% 
  log2() %>% 
  as.matrix() %>% 
  pheatmap(scale = "row",
           show_rownames = F,
           show_colnames = F,
           cluster_cols = as.hclust(col_dend),
           cluster_rows = F,
           annotation = annotation,
           col = inferno(100))
dev.off()
