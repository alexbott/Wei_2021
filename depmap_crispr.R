# Author: Alex Bott
# Function: Analyze DLBCL dependencies from DepMap
# Date: 06/17/21

# imports
library(tidyverse)
library(data.table)

# data import
meta <- read.csv("data/depmap21Q2/sample_info.csv")
crispr <- fread("data/depmap21Q2/CRISPR_gene_effect.csv")

# 
slc25a1 <- data.frame(DepMap_ID = crispr$DepMap_ID,
                      slc = crispr$`SLC25A1 (6576)`) %>% 
  left_join(meta)

# For every lineage in our data set with at least five cell lines, we calculated 
# the difference in means in gene dependency between cell lines of that lineage and 
# the rest of the data set, and assessed significance with a one-tailed 
# Student's t test (df = 340), for each gene screened. Differential dependencies 
# were called with a negative effect size at a significance of FDR-corrected P 
# value <0.05. For each chromosome arm that was recurrently amplified for that lineage, 
# we calculated the fraction of significant differential dependencies on that 
# chromosome arm before and after CERES correction.
# Meyers et al 2017

# lymphocyte lineage
slc25a1 %>% 
  filter(!(is.na(lineage))) %>% 
  filter(!(grepl("engineered", lineage))) %>% 
  mutate(lineage = gsub("_", " ", lineage)) %>% 
  mutate(lineage = str_to_title(lineage)) %>% 
  ggplot(aes(x = slc, y = lineage)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.35) + 
  scale_y_discrete(limits = rev) +
  labs(x = "Gene Effect (CERES)", y = "") +
  theme_bw()

slc25a1 %>% 
  filter(!(is.na(lineage))) %>% 
  filter(!(grepl("engineered", lineage))) %>% 
  mutate(lineage = gsub("_", " ", lineage)) %>% 
  mutate(lineage = str_to_title(lineage)) %>% 
  mutate(dlbcl = case_when(lineage_sub_subtype == "DLBCL" ~ TRUE,
                           TRUE ~ FALSE)) %>% 
  ggplot(aes(x = slc, y = lineage)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(col = dlbcl),
             alpha = 0.35) + 
  scale_y_discrete(limits = rev) +
  scale_color_manual(name = "Subtype", 
                     labels = c("Other", "DLBCL"),
                     values = c("black", "magenta")) +
  labs(x = "Gene Effect (CERES)", y = "") +
  theme_bw()

slc25a1 %>% 
  mutate(lymph = str_detect(lineage, "lymphocyte")) %>% 
  summarize(lymphocyteMean = mean(slc[lymph], na.rm = T),
            otherMean = mean(slc[!lymph], na.rm = T),
            diff = lymphocyteMean - otherMean,
            pvalue = t.test(slc[lymph], slc[!lymph],
                            var.equal = T, alternative = "two.sided")$p.value)

# DLBCL subtype
slc25a1 %>% 
  filter(!(is.na(lineage))) %>% 
  filter(!(grepl("engineered", lineage))) %>% 
  mutate(lineage = gsub("_", " ", lineage)) %>% 
  mutate(lineage = str_to_title(lineage)) %>% 
  ggplot(aes(x = slc, y = lineage_sub_subtype)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.35) + 
  scale_y_discrete(limits = rev) +
  labs(x = "Gene Effect (CERES)", y = "") +
  theme_bw()

slc25a1 %>% 
  mutate(dlbcl = str_detect(lineage_sub_subtype, "DLBCL")) %>% 
  summarize(dlbclMean = mean(slc[dlbcl], na.rm = T),
            otherMean = mean(slc[!dlbcl], na.rm = T),
            diff = dlbclMean - otherMean,
            pvalue = t.test(slc[dlbcl], slc[!dlbcl],
                            var.equal = T, alternative = "two.sided")$p.value)
