setwd("C:/Users/cnewt/Yang/Dawei_RNA/counts") #For local computer
work <- ("/work/lylab/cjn40747/dawei_RNA/counts") #For Sapelo2
setwd(work)

#This script uses the conda environment "heatmap"
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(tidyverse)
library(dplyr)
library(ggvenn)
library(DESeq2)
library(tximport)
library(readr)
library(rhdf5)
library(ensembldb)
library(plyr)
library(ballgown)
library(txdbmaker)
library(purrr)

set.seed(1)

#Data loading and prep----
samples <- read.csv('C:/Users/cnewt/Yang/Dawei_RNA/counts/h5_files_meta.csv')
rownames(samples) <- samples$sample
files <- file.path('C:/Users/cnewt/Yang/Dawei_RNA/h5_files', samples$sample, 'abundance.h5')
txdb <- makeTxDbFromGFF('C:/Users/cnewt/Yang/Dawei_RNA/counts/TAIR10_DNA.gtf')
k <- keys(txdb, keytype = 'TXNAME')
tx_map <- select(txdb, keys = k, columns = 'GENEID', keytype = 'TXNAME')
tx2gene <- tx_map
tx2gene$TXNAME <- gsub('transcript:', '', tx2gene$TXNAME)
tx2gene$GENEID <- gsub('gene:', '', tx2gene$GENEID)
txi.kallisto <- tximport(files, type = 'kallisto', tx2gene = tx2gene, ignoreAfterBar = TRUE)
meta <- read.csv('meta.csv')
str(meta)
meta$Genotype <- as.factor(meta$Genotype)
meta$Time <- as.factor(meta$Time)
Group <- factor(paste(meta$Genotype,meta$Time,sep='.'))
meta_compiled <- cbind(meta, Group=Group)
meta_compiled <- model.matrix(~0+Group)
colnames(meta_compiled) <- levels(Group)

hour0 <- as.matrix(compiled[, c(1:3, 22:24, 43:45, 64:66)])
min10 <- as.matrix(compiled[, c(4:6, 25:27, 46:48, 67:69)])
hour1 <- as.matrix(compiled[, c(7:9, 28:30, 49:51, 70:72)])
hour24 <- as.matrix(compiled[, c(10:12, 31:33, 52:54, 73:75)])
min30 <- as.matrix(compiled[, c(13:15, 34:36, 55:57, 76:78)])
hour48 <- as.matrix(compiled[, c(16:18, 37:39, 58:60, 79:81)])
hour6 <- as.matrix(compiled[, c(19:21, 40:42, 61:63, 82:84)])

#All treatments----
counts <- DGEList(txi.kallisto$counts, samples = meta)
keep_compiled <- filterByExpr(counts, meta_compiled) 
counts <- counts[keep_compiled, , keep.lib.sizes=FALSE]
counts <- calcNormFactors(counts, method='TMM')
counts <- normLibSizes(counts)
counts_matrix <- as.matrix(counts$counts)
str(counts_matrix)
levels <- c("Col.0.1", "Col.0.2", "Col.0.3", "npr3.4.0.1", "npr3.4.0.2", "npr3.4.0.3", "tga256.0.1", "tga256.0.2", "tga256.0.3", "sid2.0.1", "sid2.0.2", "sid2.0.3", "Col.10m.1", "Col.10m.2", "Col.10m.3", "npr3.4.10m.1", "npr3.4.10m.2", "npr3.4.10m.3", "tga256.10m.1", "tga256.10m.2", "tga256.10m.3", "sid2.10m.1", "sid2.10m.2", "sid2.10m.3", "Col.30m.1", "Col.30m.2", "Col.30m.3", "npr3.4.30m.1", "npr3.4.30m.2", "npr3.4.30m.3", "tga256.30m.1", "tga256.30m.2", "tga256.30m.3", "sid2.30m.1", "sid2.30m.2", "sid2.30m.3", "Col.1hr.1", "Col.1hr.2", "Col.1hr.3", "npr3.4.1hr.1", "npr3.4.1hr.2", "npr3.4.1hr.3", "tga256.1hr.1", "tga256.1hr.2", "tga256.1hr.3", "sid2.1hr.1", "sid2.1hr.2", "sid2.1hr.3", "Col.6hr.1", "Col.6hr.2", "Col.6hr.3", "npr3.4.6hr.1", "npr3.4.6hr.2", "npr3.4.6hr.3", "tga256.6hr.1", "tga256.6hr.2", "tga256.6hr.3", "sid2.6hr.1", "sid2.6hr.2", "sid2.6hr.3", "Col.24hr.1", "Col.24hr.2", "Col.24hr.3", "npr3.4.24hr.1", "npr3.4.24hr.2", "npr3.4.24hr.3", "tga256.24hr.1", "tga256.24hr.2", "tga256.24hr.3", "sid2.24hr.1", "sid2.24hr.2", "sid2.24hr.3", "Col.48hr.1", "Col.48hr.2", "Col.48hr.3", "npr3.4.48hr.1", "npr3.4.48hr.2", "npr3.4.48hr.3", "tga256.48hr.1", "tga256.48hr.2", "tga256.48hr.3", "sid2.48hr.1", "sid2.48hr.2", "sid2.48hr.3")
transcript_heatmap <- pheatmap(counts_matrix[, levels], scale = 'row', show_rownames = FALSE, color = brewer.pal(11, 'RdYlBu'), cutree_rows = 7, clustering_method = 'ward.D', cluster_cols = FALSE, gaps_col = c(12,24,36,48,60,72,84))
jpeg('heatmap.jpg', width = 12, height = 12, units = 'in', res = 300)
print(transcript_heatmap)
dev.off()

geno_labels <- c('Col-0', substitute(paste(italic('npr34'))), substitute(paste(italic('sid2'))), substitute(paste(italic('tga256'))))
time_labels <- c('0h', '10min', '30min', '1h', '6h', '24h', '48h')
geno_group <- c(15,16,17,18)
time_group <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")

jpeg('MD_Plot.jpg', width = 6, height = 6, units = 'in', res = 300)
par(mar = c(5.1,4.1,4.1,7))
plotMDS(counts, col = time_group[counts$samples$Time], pch = geno_group[counts$samples$Genotype], xlab = 'Dimension 1', ylab = 'Dimension 2', gene.selection = 'common')
legend(legend=geno_labels, pch = geno_group, xpd = TRUE, title = 'Genotype', x = 1.4655, y = 1.4265)
legend(legend =time_labels, col = time_group, pch = 15, xpd = TRUE, title = 'Time', x = 1.4655, y = 0.665)
dev.off()

#10 min data----
match <- colnames(min10)
meta_min10 <- meta[meta$ID %in% match, ]
counts_min10 <- DGEList(min10, samples = meta_min10)
keep_min10 <- filterByExpr(counts_min10)
counts_min10 <- counts_min10[keep_min10, , keep.lib.sizes=FALSE]
counts_min10 <- calcNormFactors(counts_min10, method='TMM')
counts_min10 <- normLibSizes(counts_min10)
counts_min10_matrix <- as.matrix(counts_min10$counts)
str(counts_min10_matrix)
levels_min10 <- c("Col.10m.1", "Col.10m.2", "Col.10m.3", "npr3.4.10m.1", "npr3.4.10m.2", "npr3.4.10m.3", "tga256.10m.1", "tga256.10m.2", "tga256.10m.3", "sid2.10m.1", "sid2.10m.2", "sid2.10m.3")
transcript_heatmap_min10 <- pheatmap(counts_min10_matrix[, levels_min10], scale = 'row', show_rownames = FALSE, color = brewer.pal(11, 'RdYlBu'), clustering_method = 'ward.D', cluster_cols = FALSE)
jpeg('heatmap_min10.jpg', width = 6, height = 12, units = 'in', res = 300)
print(transcript_heatmap_min10)
dev.off()

#30 min data----
counts_min30 <- DGEList(min30)
check_min30 <- min30 > 1
keep_min30 <- which(rowSums(check_min30) >= 2)
counts_min30 <- counts_min30[keep_min30,]
counts_min30 <- calcNormFactors(counts_min30, method='TMM')
counts_min30 <- normLibSizes(counts_min30)
counts_min30_matrix <- as.matrix(counts_min30$counts)
str(counts_min30_matrix)
levels_min30 <- c("Col.30m.1", "Col.30m.2", "Col.30m.3", "npr3.4.30m.1", "npr3.4.30m.2", "npr3.4.30m.3", "tga256.30m.1", "tga256.30m.2", "tga256.30m.3", "sid2.30m.1", "sid2.30m.2", "sid2.30m.3")
transcript_heatmap_min30 <- pheatmap(counts_min30_matrix[, levels_min30], scale = 'row', show_rownames = FALSE, color = brewer.pal(11, 'RdYlBu'), clustering_method = 'ward.D', cluster_cols = FALSE)
jpeg('heatmap_min30.jpg', width = 6, height = 12, units = 'in', res = 300)
print(transcript_heatmap_min30)
dev.off()

#1 hour data----
counts_hour1 <- DGEList(hour1)
check_hour1 <- hour1 > 1
keep_hour1 <- which(rowSums(check_hour1) >= 2)
counts_hour1 <- counts_hour1[keep_hour1,]
counts_hour1 <- calcNormFactors(counts_hour1, method='TMM')
counts_hour1 <- normLibSizes(counts_hour1)
counts_hour1_matrix <- as.matrix(counts_hour1$counts)
str(counts_hour1_matrix)
levels_hour1 <- c("Col.1hr.1", "Col.1hr.2", "Col.1hr.3", "npr3.4.1hr.1", "npr3.4.1hr.2", "npr3.4.1hr.3", "tga256.1hr.1", "tga256.1hr.2", "tga256.1hr.3", "sid2.1hr.1", "sid2.1hr.2", "sid2.1hr.3")
transcript_heatmap_hour1 <- pheatmap(counts_hour1_matrix[, levels_hour1], scale = 'row', show_rownames = FALSE, color = brewer.pal(11, 'RdYlBu'), clustering_method = 'ward.D', cluster_cols = FALSE)
jpeg('heatmap_hour1.jpg', width = 6, height = 12, units = 'in', res = 300)
print(transcript_heatmap_hour1)
dev.off()

#6 hour data----
counts_hour6 <- DGEList(hour6)
check_hour6 <- hour6 > 1
keep_hour6 <- which(rowSums(check_hour6) >= 2)
counts_hour6 <- counts_hour6[keep_hour6,]
counts_hour6 <- calcNormFactors(counts_hour6, method='TMM')
counts_hour6 <- normLibSizes(counts_hour6)
counts_hour6_matrix <- as.matrix(counts_hour6$counts)
str(counts_hour6_matrix)
levels_hour6 <- c("Col.6hr.1", "Col.6hr.2", "Col.6hr.3", "npr3.4.6hr.1", "npr3.4.6hr.2", "npr3.4.6hr.3", "tga256.6hr.1", "tga256.6hr.2", "tga256.6hr.3", "sid2.6hr.1", "sid2.6hr.2", "sid2.6hr.3")
transcript_heatmap_hour6 <- pheatmap(counts_hour6_matrix[, levels_hour6], scale = 'row', show_rownames = FALSE, color = brewer.pal(11, 'RdYlBu'), clustering_method = 'ward.D', cluster_cols = FALSE)
jpeg('heatmap_hour6.jpg', width = 6, height = 12, units = 'in', res = 300)
print(transcript_heatmap_hour6)
dev.off()

#24 hour data----
counts_hour24 <- DGEList(hour24)
check_hour24 <- hour24 > 1
keep_hour24 <- which(rowSums(check_hour24) >= 2)
counts_hour24 <- counts_hour24[keep_hour24,]
counts_hour24 <- calcNormFactors(counts_hour24, method='TMM')
counts_hour24 <- normLibSizes(counts_hour24)
counts_hour24_matrix <- as.matrix(counts_hour24$counts)
str(counts_hour24_matrix)
levels_hour24 <- c("Col.24hr.1", "Col.24hr.2", "Col.24hr.3", "npr3.4.24hr.1", "npr3.4.24hr.2", "npr3.4.24hr.3", "tga256.24hr.1", "tga256.24hr.2", "tga256.24hr.3", "sid2.24hr.1", "sid2.24hr.2", "sid2.24hr.3")
transcript_heatmap_hour24 <- pheatmap(counts_hour24_matrix[, levels_hour24], scale = 'row', show_rownames = FALSE, color = brewer.pal(11, 'RdYlBu'), clustering_method = 'ward.D', cluster_cols = FALSE)
jpeg('heatmap_hour24.jpg', width = 6, height = 12, units = 'in', res = 300)
print(transcript_heatmap_hour24)
dev.off()

#48 hour data----
counts_hour48 <- DGEList(hour48)
check_hour48 <- hour48 > 1
keep_hour48 <- which(rowSums(check_hour48) >= 2)
counts_hour48 <- counts_hour48[keep_hour48,]
counts_hour48 <- calcNormFactors(counts_hour48, method='TMM')
counts_hour48 <- normLibSizes(counts_hour48)
counts_hour48_matrix <- as.matrix(counts_hour48$counts)
str(counts_hour48_matrix)
levels_hour48 <- c("Col.48hr.1", "Col.48hr.2", "Col.48hr.3", "npr3.4.48hr.1", "npr3.4.48hr.2", "npr3.4.48hr.3", "tga256.48hr.1", "tga256.48hr.2", "tga256.48hr.3", "sid2.48hr.1", "sid2.48hr.2", "sid2.48hr.3")
transcript_heatmap_hour48 <- pheatmap(counts_hour48_matrix[, levels_hour48], scale = 'row', show_rownames = FALSE, color = brewer.pal(11, 'RdYlBu'), clustering_method = 'ward.D', cluster_cols = FALSE)
jpeg('heatmap_hour48.jpg', width = 6, height = 12, units = 'in', res = 300)
print(transcript_heatmap_hour48)
dev.off()

#DEG----
counts <- estimateDisp(counts, meta_compiled)
fit_compiled <- glmQLFit(counts, meta_compiled)

compiled_contrasts <- makeContrasts(
  col10 = A.2-A.1,
  col30 = A.3-A.1,
  col1 = A.4-A.1,
  col6 = A.5-A.1,
  col24 = A.6-A.1,
  col48 = A.7-A.1,
  npr10 = B.2-B.1,
  npr30 = B.3-B.1,
  npr1 = B.4-B.1,
  npr6 = B.5-B.1,
  npr24 = B.6-B.1,
  npr48 = B.7-B.1,
  sid10 = C.2-C.1,
  sid30 = C.3-C.1,
  sid1 = C.4-C.1,
  sid6 = C.5-C.1,
  sid24 = C.6-C.1,
  sid48 = C.7-C.1,
  tga10 = D.2-D.1,
  tga30 = D.3-D.1,
  tga1 = D.4-D.1,
  tga6 = D.5-D.1,
  tga24 = D.6-D.1,
  tga48 = D.7-D.1,
  col_npr0 = B.1-A.1,
  col_npr10 = B.2-A.2,
  col_npr30 = B.3-A.3,
  col_npr1 = B.4-A.4,
  col_npr6 = B.5-A.5,
  col_npr24 = B.6-A.6,
  col_npr48 = B.7-A.7,
  col_sid0 = C.1-A.1,
  col_sid10 = C.2-A.2,
  col_sid30 = C.3-A.3,
  col_sid1 = C.4-A.4,
  col_sid6 = C.5-A.5,
  col_sid24 = C.6-A.6,
  col_sid48 = C.7-A.7,
  col_tga0 = D.1-A.1,
  col_tga10 = D.2-A.2,
  col_tga30 = D.3-A.3,
  col_tga1 = D.4-A.4,
  col_tga6 = D.5-A.5,
  col_tga24 = D.6-A.6,
  col_tga48 = D.7-A.7,
  sid_tga0 = C.1-D.1,
  sid_tga10 = C.2-D.2,
  sid_tga30 = C.3-D.3,
  sid_tga1 = C.4-D.4,
  sid_tga6 = C.5-D.5,
  sid_tga24 = C.6-D.6,
  sid_tga48 = C.7-D.7,
  levels = meta_compiled)

##Analysis 1----
###24 hour----
A1_col24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col24'])
A1_col24_deg_list <- topTags(A1_col24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_col24_deg_df <- A1_col24_deg_list$table
A1_col24_deg_df <- A1_col24_deg_df[A1_col24_deg_df$logFC > 0.5, ]
A1_col24_deg_df <- rownames_to_column(A1_col24_deg_df, var = 'Gene')
names(A1_col24_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A1_sid24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid24'])
A1_sid24_deg_list <- topTags(A1_sid24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_sid24_deg_df <- A1_sid24_deg_list$table
A1_sid24_deg_df <- A1_sid24_deg_df[A1_sid24_deg_df$logFC > 0.5, ]
A1_sid24_deg_df <- rownames_to_column(A1_sid24_deg_df, var = 'Gene')
names(A1_sid24_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A1_tga24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga24'])
A1_tga24_deg_list <- topTags(A1_tga24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_tga24_deg_df <- A1_tga24_deg_list$table
A1_tga24_deg_df <- A1_tga24_deg_df[A1_tga24_deg_df$logFC < -0.5, ]
A1_tga24_deg_df <- rownames_to_column(A1_tga24_deg_df, var = 'Gene')
names(A1_tga24_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A1_DEG_list_24 <- list(A1_col24_deg_df, A1_sid24_deg_df, A1_tga24_deg_df)
A1_DEG_compiled_24 <- A1_DEG_list_24 %>% purrr::reduce(full_join, by = 'Gene')
A1_DEG_compiled_24 <- na.omit(A1_DEG_compiled_24)
write.csv(A1_DEG_compiled_24, 'A1_24_deg_df.csv')


###10 min----
A1_col10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col10'])
A1_col10_deg_list <- topTags(A1_col10_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_col10_deg_df <- A1_col10_deg_list$table
A1_col10_deg_df <- A1_col10_deg_df[A1_col10_deg_df$logFC > 0.5, ]
A1_col10_deg_df <- rownames_to_column(A1_col10_deg_df, var = 'Gene')
names(A1_col10_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A1_sid10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid10'])
A1_sid10_deg_list <- topTags(A1_sid10_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_sid10_deg_df <- A1_sid10_deg_list$table
A1_sid10_deg_df <- A1_sid10_deg_df[A1_sid10_deg_df$logFC > 0.5, ]
A1_sid10_deg_df <- rownames_to_column(A1_sid10_deg_df, var = 'Gene')
names(A1_sid10_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A1_tga10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga10'])
A1_tga10_deg_list <- topTags(A1_tga10_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_tga10_deg_df <- A1_tga10_deg_list$table
A1_tga10_deg_df <- A1_tga10_deg_df[A1_tga10_deg_df$logFC < -0.5, ]
A1_tga10_deg_df <- rownames_to_column(A1_tga10_deg_df, var = 'Gene')
names(A1_tga10_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A1_DEG_list_10 <- list(A1_col10_deg_df, A1_sid10_deg_df, A1_tga10_deg_df)
A1_DEG_compiled_10 <- A1_DEG_list_10 %>% purrr::reduce(full_join, by = 'Gene')
A1_DEG_compiled_10 <- na.omit(A1_DEG_compiled_10)
write.csv(A1_DEG_compiled_10, 'A1_10_deg_df.csv')

###30 min----
A1_col30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col30'])
A1_col30_deg_list <- topTags(A1_col30_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_col30_deg_df <- A1_col30_deg_list$table
A1_col30_deg_df <- A1_col30_deg_df[A1_col30_deg_df$logFC > 0.5, ]
A1_col30_deg_df <- rownames_to_column(A1_col30_deg_df, var = 'Gene')
names(A1_col30_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A1_sid30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid30'])
A1_sid30_deg_list <- topTags(A1_sid30_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_sid30_deg_df <- A1_sid30_deg_list$table
A1_sid30_deg_df <- A1_sid30_deg_df[A1_sid30_deg_df$logFC > 0.5, ]
A1_sid30_deg_df <- rownames_to_column(A1_sid30_deg_df, var = 'Gene')
names(A1_sid30_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A1_tga30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga30'])
A1_tga30_deg_list <- topTags(A1_tga30_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_tga30_deg_df <- A1_tga30_deg_list$table
A1_tga30_deg_df <- A1_tga30_deg_df[A1_tga30_deg_df$logFC < -0.5, ]
A1_tga30_deg_df <- rownames_to_column(A1_tga30_deg_df, var = 'Gene')
names(A1_tga30_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A1_DEG_list_30 <- list(A1_col30_deg_df, A1_sid30_deg_df, A1_tga30_deg_df)
A1_DEG_compiled_30 <- A1_DEG_list_30 %>% purrr::reduce(full_join, by = 'Gene')
A1_DEG_compiled_30 <- na.omit(A1_DEG_compiled_30)
write.csv(A1_DEG_compiled_30, 'A1_30_deg_df.csv')

###1 hour----
A1_col1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col1'])
A1_col1_deg_list <- topTags(A1_col1_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_col1_deg_df <- A1_col1_deg_list$table
A1_col1_deg_df <- A1_col1_deg_df[A1_col1_deg_df$logFC > 0.5, ]
A1_col1_deg_df <- rownames_to_column(A1_col1_deg_df, var = 'Gene')
names(A1_col1_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A1_sid1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid1'])
A1_sid1_deg_list <- topTags(A1_sid1_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_sid1_deg_df <- A1_sid1_deg_list$table
A1_sid1_deg_df <- A1_sid1_deg_df[A1_sid1_deg_df$logFC > 0.5, ]
A1_sid1_deg_df <- rownames_to_column(A1_sid1_deg_df, var = 'Gene')
names(A1_sid1_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A1_tga1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga1'])
A1_tga1_deg_list <- topTags(A1_tga1_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_tga1_deg_df <- A1_tga1_deg_list$table
A1_tga1_deg_df <- A1_tga1_deg_df[A1_tga1_deg_df$logFC < -0.5, ]
A1_tga1_deg_df <- rownames_to_column(A1_tga1_deg_df, var = 'Gene')
names(A1_tga1_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A1_DEG_list_1 <- list(A1_col1_deg_df, A1_sid1_deg_df, A1_tga1_deg_df)
A1_DEG_compiled_1 <- A1_DEG_list_1 %>% purrr::reduce(full_join, by = 'Gene')
A1_DEG_compiled_1 <- na.omit(A1_DEG_compiled_1)
write.csv(A1_DEG_compiled_1, 'A1_1_deg_df.csv')

###6 hour----
A1_col6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col6'])
A1_col6_deg_list <- topTags(A1_col6_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_col6_deg_df <- A1_col6_deg_list$table
A1_col6_deg_df <- A1_col6_deg_df[A1_col6_deg_df$logFC > 0.5, ]
A1_col6_deg_df <- rownames_to_column(A1_col6_deg_df, var = 'Gene')
names(A1_col6_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A1_sid6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid6'])
A1_sid6_deg_list <- topTags(A1_sid6_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_sid6_deg_df <- A1_sid6_deg_list$table
A1_sid6_deg_df <- A1_sid6_deg_df[A1_sid6_deg_df$logFC > 0.5, ]
A1_sid6_deg_df <- rownames_to_column(A1_sid6_deg_df, var = 'Gene')
names(A1_sid6_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A1_tga6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga6'])
A1_tga6_deg_list <- topTags(A1_tga6_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_tga6_deg_df <- A1_tga6_deg_list$table
A1_tga6_deg_df <- A1_tga6_deg_df[A1_tga6_deg_df$logFC < -0.5, ]
A1_tga6_deg_df <- rownames_to_column(A1_tga6_deg_df, var = 'Gene')
names(A1_tga6_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A1_DEG_list_6 <- list(A1_col6_deg_df, A1_sid6_deg_df, A1_tga6_deg_df)
A1_DEG_compiled_6 <- A1_DEG_list_6 %>% purrr::reduce(full_join, by = 'Gene')
A1_DEG_compiled_6 <- na.omit(A1_DEG_compiled_6)
write.csv(A1_DEG_compiled_6, 'A1_6_deg_df.csv')

###48 hour----
A1_col48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col48'])
A1_col48_deg_list <- topTags(A1_col48_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_col48_deg_df <- A1_col48_deg_list$table
A1_col48_deg_df <- A1_col48_deg_df[A1_col48_deg_df$logFC > 0.5, ]
A1_col48_deg_df <- rownames_to_column(A1_col48_deg_df, var = 'Gene')
names(A1_col48_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A1_sid48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid48'])
A1_sid48_deg_list <- topTags(A1_sid48_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_sid48_deg_df <- A1_sid48_deg_list$table
A1_sid48_deg_df <- A1_sid48_deg_df[A1_sid48_deg_df$logFC > 0.5, ]
A1_sid48_deg_df <- rownames_to_column(A1_sid48_deg_df, var = 'Gene')
names(A1_sid48_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A1_tga48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga48'])
A1_tga48_deg_list <- topTags(A1_tga48_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A1_tga48_deg_df <- A1_tga48_deg_list$table
A1_tga48_deg_df <- A1_tga48_deg_df[A1_tga48_deg_df$logFC < -0.5, ]
A1_tga48_deg_df <- rownames_to_column(A1_tga48_deg_df, var = 'Gene')
names(A1_tga48_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A1_DEG_list_48 <- list(A1_col48_deg_df, A1_sid48_deg_df, A1_tga48_deg_df)
A1_DEG_compiled_48 <- A1_DEG_list_48 %>% purrr::reduce(full_join, by = 'Gene')
A1_DEG_compiled_48 <- na.omit(A1_DEG_compiled_48)
write.csv(A1_DEG_compiled_48, 'A1_48_deg_df.csv')

##Analysis 2----
###24 hour----
A2_col24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col24'])
A2_col24_deg_list <- topTags(A2_col24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_col24_deg_df <- A2_col24_deg_list$table
A2_col24_deg_df <- A2_col24_deg_df[A2_col24_deg_df$logFC < -0.5, ]
A2_col24_deg_df <- rownames_to_column(A2_col24_deg_df, var = 'Gene')
names(A2_col24_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A2_sid24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid24'])
A2_sid24_deg_list <- topTags(A2_sid24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_sid24_deg_df <- A2_sid24_deg_list$table
A2_sid24_deg_df <- A2_sid24_deg_df[A2_sid24_deg_df$logFC < -0.5, ]
A2_sid24_deg_df <- rownames_to_column(A2_sid24_deg_df, var = 'Gene')
names(A2_sid24_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A2_tga24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga24'])
A2_tga24_deg_list <- topTags(A2_tga24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_tga24_deg_df <- A2_tga24_deg_list$table
A2_tga24_deg_df <- A2_tga24_deg_df[A2_tga24_deg_df$logFC > 0.5, ]
A2_tga24_deg_df <- rownames_to_column(A2_tga24_deg_df, var = 'Gene')
names(A2_tga24_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A2_DEG_list_24 <- list(A2_col24_deg_df, A2_sid24_deg_df, A2_tga24_deg_df)
A2_DEG_compiled_24 <- A2_DEG_list_24 %>% purrr::reduce(full_join, by = 'Gene')
A2_DEG_compiled_24 <- na.omit(A2_DEG_compiled_24)
write.csv(A2_DEG_compiled_24, 'A2_24_deg_df.csv')

###10 min----
A2_col10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col10'])
A2_col10_deg_list <- topTags(A2_col10_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_col10_deg_df <- A2_col10_deg_list$table
A2_col10_deg_df <- A2_col10_deg_df[A2_col10_deg_df$logFC < -0.5, ]
A2_col10_deg_df <- rownames_to_column(A2_col10_deg_df, var = 'Gene')
names(A2_col10_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A2_sid10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid10'])
A2_sid10_deg_list <- topTags(A2_sid10_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_sid10_deg_df <- A2_sid10_deg_list$table
A2_sid10_deg_df <- A2_sid10_deg_df[A2_sid10_deg_df$logFC < -0.5, ]
A2_sid10_deg_df <- rownames_to_column(A2_sid10_deg_df, var = 'Gene')
names(A2_sid10_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A2_tga10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga10'])
A2_tga10_deg_list <- topTags(A2_tga10_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_tga10_deg_df <- A2_tga10_deg_list$table
A2_tga10_deg_df <- A2_tga10_deg_df[A2_tga10_deg_df$logFC > 0.5, ]
A2_tga10_deg_df <- rownames_to_column(A2_tga10_deg_df, var = 'Gene')
names(A2_tga10_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A2_DEG_list_10 <- list(A2_col10_deg_df, A2_sid10_deg_df, A2_tga10_deg_df)
A2_DEG_compiled_10 <- A2_DEG_list_10 %>% purrr::reduce(full_join, by = 'Gene')
A2_DEG_compiled_10 <- na.omit(A2_DEG_compiled_10)
write.csv(A2_DEG_compiled_10, 'A2_10_deg_df.csv')

###30 min----
A2_col30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col30'])
A2_col30_deg_list <- topTags(A2_col30_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_col30_deg_df <- A2_col30_deg_list$table
A2_col30_deg_df <- A2_col30_deg_df[A2_col30_deg_df$logFC < -0.5, ]
A2_col30_deg_df <- rownames_to_column(A2_col30_deg_df, var = 'Gene')
names(A2_col30_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A2_sid30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid30'])
A2_sid30_deg_list <- topTags(A2_sid30_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_sid30_deg_df <- A2_sid30_deg_list$table
A2_sid30_deg_df <- A2_sid30_deg_df[A2_sid30_deg_df$logFC < -0.5, ]
A2_sid30_deg_df <- rownames_to_column(A2_sid30_deg_df, var = 'Gene')
names(A2_sid30_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A2_tga30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga30'])
A2_tga30_deg_list <- topTags(A2_tga30_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_tga30_deg_df <- A2_tga30_deg_list$table
A2_tga30_deg_df <- A2_tga30_deg_df[A2_tga30_deg_df$logFC > 0.5, ]
A2_tga30_deg_df <- rownames_to_column(A2_tga30_deg_df, var = 'Gene')
names(A2_tga30_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A2_DEG_list_30 <- list(A2_col30_deg_df, A2_sid30_deg_df, A2_tga30_deg_df)
A2_DEG_compiled_30 <- A2_DEG_list_30 %>% purrr::reduce(full_join, by = 'Gene')
A2_DEG_compiled_30 <- na.omit(A2_DEG_compiled_30)
write.csv(A2_DEG_compiled_30, 'A2_30_deg_df.csv')

###1 hour----
A2_col1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col1'])
A2_col1_deg_list <- topTags(A2_col1_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_col1_deg_df <- A2_col1_deg_list$table
A2_col1_deg_df <- A2_col1_deg_df[A2_col1_deg_df$logFC < -0.5, ]
A2_col1_deg_df <- rownames_to_column(A2_col1_deg_df, var = 'Gene')
names(A2_col1_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A2_sid1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid1'])
A2_sid1_deg_list <- topTags(A2_sid1_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_sid1_deg_df <- A2_sid1_deg_list$table
A2_sid1_deg_df <- A2_sid1_deg_df[A2_sid1_deg_df$logFC < -0.5, ]
A2_sid1_deg_df <- rownames_to_column(A2_sid1_deg_df, var = 'Gene')
names(A2_sid1_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A2_tga1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga1'])
A2_tga1_deg_list <- topTags(A2_tga1_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_tga1_deg_df <- A2_tga1_deg_list$table
A2_tga1_deg_df <- A2_tga1_deg_df[A2_tga1_deg_df$logFC > 0.5, ]
A2_tga1_deg_df <- rownames_to_column(A2_tga1_deg_df, var = 'Gene')
names(A2_tga1_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A2_DEG_list_1 <- list(A2_col1_deg_df, A2_sid1_deg_df, A2_tga1_deg_df)
A2_DEG_compiled_1 <- A2_DEG_list_1 %>% purrr::reduce(full_join, by = 'Gene')
A2_DEG_compiled_1 <- na.omit(A2_DEG_compiled_1)
write.csv(A2_DEG_compiled_1, 'A2_1_deg_df.csv')

###6 hour----
A2_col6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col6'])
A2_col6_deg_list <- topTags(A2_col6_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_col6_deg_df <- A2_col6_deg_list$table
A2_col6_deg_df <- A2_col6_deg_df[A2_col6_deg_df$logFC < -0.5, ]
A2_col6_deg_df <- rownames_to_column(A2_col6_deg_df, var = 'Gene')
names(A2_col6_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A2_sid6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid6'])
A2_sid6_deg_list <- topTags(A2_sid6_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_sid6_deg_df <- A2_sid6_deg_list$table
A2_sid6_deg_df <- A2_sid6_deg_df[A2_sid6_deg_df$logFC < -0.5, ]
A2_sid6_deg_df <- rownames_to_column(A2_sid6_deg_df, var = 'Gene')
names(A2_sid6_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A2_tga6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga6'])
A2_tga6_deg_list <- topTags(A2_tga6_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_tga6_deg_df <- A2_tga6_deg_list$table
A2_tga6_deg_df <- A2_tga6_deg_df[A2_tga6_deg_df$logFC > 0.5, ]
A2_tga6_deg_df <- rownames_to_column(A2_tga6_deg_df, var = 'Gene')
names(A2_tga6_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A2_DEG_list_6 <- list(A2_col6_deg_df, A2_sid6_deg_df, A2_tga6_deg_df)
A2_DEG_compiled_6 <- A2_DEG_list_6 %>% purrr::reduce(full_join, by = 'Gene')
A2_DEG_compiled_6 <- na.omit(A2_DEG_compiled_6)
write.csv(A2_DEG_compiled_6, 'A2_6_deg_df.csv')

###48 hour----
A2_col48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col48'])
A2_col48_deg_list <- topTags(A2_col48_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_col48_deg_df <- A2_col48_deg_list$table
A2_col48_deg_df <- A2_col48_deg_df[A2_col48_deg_df$logFC < -0.5, ]
A2_col48_deg_df <- rownames_to_column(A2_col48_deg_df, var = 'Gene')
names(A2_col48_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')

A2_sid48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid48'])
A2_sid48_deg_list <- topTags(A2_sid48_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_sid48_deg_df <- A2_sid48_deg_list$table
A2_sid48_deg_df <- A2_sid48_deg_df[A2_sid48_deg_df$logFC < -0.5, ]
A2_sid48_deg_df <- rownames_to_column(A2_sid48_deg_df, var = 'Gene')
names(A2_sid48_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')

A2_tga48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga48'])
A2_tga48_deg_list <- topTags(A2_tga48_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A2_tga48_deg_df <- A2_tga48_deg_list$table
A2_tga48_deg_df <- A2_tga48_deg_df[A2_tga48_deg_df$logFC > 0.5, ]
A2_tga48_deg_df <- rownames_to_column(A2_tga48_deg_df, var = 'Gene')
names(A2_tga48_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')

A2_DEG_list_48 <- list(A2_col48_deg_df, A2_sid48_deg_df, A2_tga48_deg_df)
A2_DEG_compiled_48 <- A2_DEG_list_48 %>% purrr::reduce(full_join, by = 'Gene')
A2_DEG_compiled_48 <- na.omit(A2_DEG_compiled_48)
write.csv(A2_DEG_compiled_48, 'A2_48_deg_df.csv')

##Analysis 3----
###0 hour----
A3_0_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[, 'sid_tga0'])
A3_0_deg_list <- topTags(A3_0_deg, sort.by = 'logFC', n= 'Inf', p.value = 0.05)
A3_0_deg_df <- A3_0_deg_list$table
A3_0_deg_df <- rownames_to_column(A3_0_deg_df, var = 'Gene')
A3_0_deg_df_up <- A3_0_deg_df[A3_0_deg_df$logFC > 1, ]
write.csv(A3_0_deg_df_up, 'A3_0_deg_df_up.csv')
A3_0_deg_df_down <- A3_0_deg_df[A3_0_deg_df$logFC < -1, ]
write.csv(A3_0_deg_df_down, 'A3_0_deg_df_down.csv')
A3_0_deg_df_up_GO <- read.csv('A3_0_deg_df_up_GO.csv')
str(A3_0_deg_df_up_GO)

A3_0_deg_df_up_GO_fig <- ggplot(A3_0_deg_df_up_GO, aes(x = Fold_Enrichment, y = reorder(GO_BP, Fold_Enrichment))) + geom_point(mapping = aes(size = Number_Sample, color = FDR)) + theme_classic()
A3_0_deg_df_up_GO_fig

###30 min----
A3_30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[, 'sid_tga30'])
A3_30_deg_list <- topTags(A3_30_deg, sort.by = 'logFC', n= 'Inf', p.value = 0.05)
A3_30_deg_df <- A3_30_deg_list$table
A3_30_deg_df <- rownames_to_column(A3_30_deg_df, var = 'Gene')
A3_30_deg_df_up <- A3_30_deg_df[A3_30_deg_df$logFC > 1, ]
write.csv(A3_30_deg_df_up, 'A3_30_deg_df_up.csv')
A3_30_deg_df_down <- A3_30_deg_df[A3_30_deg_df$logFC < -1, ]
write.csv(A3_30_deg_df_down, 'A3_30_deg_df_down.csv')

###1 hour----
A3_1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[, 'sid_tga1'])
A3_1_deg_list <- topTags(A3_1_deg, sort.by = 'logFC', n= 'Inf', p.value = 0.05)
A3_1_deg_df <- A3_1_deg_list$table
A3_1_deg_df <- rownames_to_column(A3_1_deg_df, var = 'Gene')
A3_1_deg_df_up <- A3_1_deg_df[A3_1_deg_df$logFC > 1, ]
write.csv(A3_1_deg_df_up, 'A3_1_deg_df_up.csv')
A3_1_deg_df_down <- A3_1_deg_df[A3_1_deg_df$logFC < -1, ]
write.csv(A3_1_deg_df_down, 'A3_1_deg_df_down.csv')

###6 hour----
A3_6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[, 'sid_tga6'])
A3_6_deg_list <- topTags(A3_6_deg, sort.by = 'logFC', n= 'Inf', p.value = 0.05)
A3_6_deg_df <- A3_6_deg_list$table
A3_6_deg_df <- rownames_to_column(A3_6_deg_df, var = 'Gene')
A3_6_deg_df_up <- A3_6_deg_df[A3_6_deg_df$logFC > 1, ]
write.csv(A3_6_deg_df_up, 'A3_6_deg_df_up.csv')
A3_6_deg_df_down <- A3_6_deg_df[A3_6_deg_df$logFC < -1, ]
write.csv(A3_6_deg_df_down, 'A3_6_deg_df_down.csv')

###24 hour----
A3_24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[, 'sid_tga24'])
A3_24_deg_list <- topTags(A3_24_deg, sort.by = 'logFC', n= 'Inf', p.value = 0.05)
A3_24_deg_df <- A3_24_deg_list$table
A3_24_deg_df <- rownames_to_column(A3_24_deg_df, var = 'Gene')
A3_24_deg_df_up <- A3_24_deg_df[A3_24_deg_df$logFC > 1, ]
write.csv(A3_24_deg_df_up, 'A3_24_deg_df_up.csv')
A3_24_deg_df_down <- A3_24_deg_df[A3_24_deg_df$logFC < -1, ]
write.csv(A3_24_deg_df_down, 'A3_24_deg_df_down.csv')

###48 hour----
A3_48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[, 'sid_tga48'])
A3_48_deg_list <- topTags(A3_48_deg, sort.by = 'logFC', n= 'Inf', p.value = 0.05)
A3_48_deg_df <- A3_48_deg_list$table
A3_48_deg_df <- rownames_to_column(A3_48_deg_df, var = 'Gene')
A3_48_deg_df_up <- A3_48_deg_df[A3_48_deg_df$logFC > 1, ]
write.csv(A3_48_deg_df_up, 'A3_48_deg_df_up.csv')
A3_48_deg_df_down <- A3_48_deg_df[A3_48_deg_df$logFC < -1, ]
write.csv(A3_48_deg_df_down, 'A3_48_deg_df_down.csv')

##Analysis 4----
###10 min----
A4_10_deg_col_df <- A1_col10_deg_df[A1_col10_deg_df$logFC > 1 | A1_col10_deg_df$logFC < -1, ]
A4_10_deg_sid_df <- A1_sid10_deg_df[A1_sid10_deg_df$logFC > 1 | A1_sid10_deg_df$logFC < -1, ]
A4_10_deg_tga_df <- A1_tga10_deg_df[A1_tga10_deg_df$logFC > 1 | A1_tga10_deg_df$logFC < -1, ]
A4_10_deg_colsid_df <- merge(A4_10_deg_col_df, A4_10_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_10_deg_colsid_notga_df <- anti_join(A4_10_deg_colsid_df, A4_10_deg_tga_df, by = 'Gene')
names(A4_10_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_10_deg_colsid_notga_df, 'A4_10_deg_df.csv')

###30 min----
A4_30_deg_col_df <- A1_col30_deg_df[A1_col30_deg_df$logFC > 1 | A1_col30_deg_df$logFC < -1, ]
A4_30_deg_sid_df <- A1_sid30_deg_df[A1_sid30_deg_df$logFC > 1 | A1_sid30_deg_df$logFC < -1, ]
A4_30_deg_tga_df <- A1_tga30_deg_df[A1_tga30_deg_df$logFC > 1 | A1_tga30_deg_df$logFC < -1, ]
A4_30_deg_colsid_df <- merge(A4_30_deg_col_df, A4_30_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_30_deg_colsid_notga_df <- anti_join(A4_30_deg_colsid_df, A4_30_deg_tga_df, by = 'Gene')
names(A4_30_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_30_deg_colsid_notga_df, 'A4_30_deg_df.csv')

###1 hour----
A4_1_deg_col_df <- A1_col1_deg_df[A1_col1_deg_df$logFC > 1 | A1_col1_deg_df$logFC < -1, ]
A4_1_deg_sid_df <- A1_sid1_deg_df[A1_sid1_deg_df$logFC > 1 | A1_sid1_deg_df$logFC < -1, ]
A4_1_deg_tga_df <- A1_tga1_deg_df[A1_tga1_deg_df$logFC > 1 | A1_tga1_deg_df$logFC < -1, ]
A4_1_deg_colsid_df <- merge(A4_1_deg_col_df, A4_1_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_1_deg_colsid_notga_df <- anti_join(A4_1_deg_colsid_df, A4_1_deg_tga_df, by = 'Gene')
names(A4_1_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_1_deg_colsid_notga_df, 'A4_1_deg_df.csv')

###6 hour----
A4_6_deg_col_df <- A1_col6_deg_df[A1_col6_deg_df$logFC > 1 | A1_col6_deg_df$logFC < -1, ]
A4_6_deg_sid_df <- A1_sid6_deg_df[A1_sid6_deg_df$logFC > 1 | A1_sid6_deg_df$logFC < -1, ]
A4_6_deg_tga_df <- A1_tga6_deg_df[A1_tga6_deg_df$logFC > 1 | A1_tga6_deg_df$logFC < -1, ]
A4_6_deg_colsid_df <- merge(A4_6_deg_col_df, A4_6_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_6_deg_colsid_notga_df <- anti_join(A4_6_deg_colsid_df, A4_6_deg_tga_df, by = 'Gene')
names(A4_6_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_6_deg_colsid_notga_df, 'A4_6_deg_df.csv')

###24 hour----
A4_24_deg_col_df <- A1_col24_deg_df[A1_col24_deg_df$logFC > 1 | A1_col24_deg_df$logFC < -1, ]
A4_24_deg_sid_df <- A1_sid24_deg_df[A1_sid24_deg_df$logFC > 1 | A1_sid24_deg_df$logFC < -1, ]
A4_24_deg_tga_df <- A1_tga24_deg_df[A1_tga24_deg_df$logFC > 1 | A1_tga24_deg_df$logFC < -1, ]
A4_24_deg_colsid_df <- merge(A4_24_deg_col_df, A4_24_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_24_deg_colsid_notga_df <- anti_join(A4_24_deg_colsid_df, A4_24_deg_tga_df, by = 'Gene')
names(A4_24_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_24_deg_colsid_notga_df, 'A4_24_deg_df.csv')

###48 hour----
A4_48_deg_col_df <- A1_col48_deg_df[A1_col48_deg_df$logFC > 1 | A1_col48_deg_df$logFC < -1, ]
A4_48_deg_sid_df <- A1_sid48_deg_df[A1_sid48_deg_df$logFC > 1 | A1_sid48_deg_df$logFC < -1, ]
A4_48_deg_tga_df <- A1_tga48_deg_df[A1_tga48_deg_df$logFC > 1 | A1_tga48_deg_df$logFC < -1, ]
A4_48_deg_colsid_df <- merge(A4_48_deg_col_df, A4_48_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_48_deg_colsid_notga_df <- anti_join(A4_48_deg_colsid_df, A4_48_deg_tga_df, by = 'Gene')
names(A4_48_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_48_deg_colsid_notga_df, 'A4_48_deg_df.csv')

##Analysis 5: Misc Comparisons----
A5_col24_deg_df <- A1_col24_deg_df[A1_col24_deg_df$logFC > 1 | A1_col24_deg_df$logFC < -1, ]
A5_col24_deg_df <- rownames_to_column(A5_col24_deg_df, var = 'Gene')
write.csv(A5_col24_deg_df, 'A5_col24_deg_df.csv')

A5_col_tga24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col_tga24'])
A5_col_tga24_deg_list <- topTags(A5_col_tga24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A5_col_tga24_deg_df <- A5_col_tga24_deg_list$table
A5_col_tga24_deg_df <- A5_col_tga24_deg_df[A5_col_tga24_deg_df$logFC > 1 | A5_col_tga24_deg_df$logFC < -1, ]
A5_col_tga24_deg_df <- rownames_to_column(A5_col_tga24_deg_df, var = 'Gene')
write.csv(A5_col_tga24_deg_df, 'A5_col_tga24_deg_df.csv')

A5_col_sid24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col_sid24'])
A5_col_sid24_deg_list <- topTags(A5_col_sid24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A5_col_sid24_deg_df <- A5_col_sid24_deg_list$table
A5_col_sid24_deg_df <- A5_col_sid24_deg_df[A5_col_sid24_deg_df$logFC > 1 | A5_col_sid24_deg_df$logFC < -1, ]
A5_col_sid24_deg_df <- rownames_to_column(A5_col_sid24_deg_df, var = 'Gene')
write.csv(A5_col_sid24_deg_df, 'A5_col_sid24_deg_df.csv')

A5_sid_tga24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid_tga24'])
A5_sid_tga24_deg_list <- topTags(A5_sid_tga24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
A5_sid_tga24_deg_df <- A5_sid_tga24_deg_list$table
A5_sid_tga24_deg_df <- A5_sid_tga24_deg_df[A5_sid_tga24_deg_df$logFC > 1 | A5_sid_tga24_deg_df$logFC < -1, ]
A5_sid_tga24_deg_df <- rownames_to_column(A5_sid_tga24_deg_df, var = 'Gene')
write.csv(A5_sid_tga24_deg_df, 'A5_sid_tga24_deg_df.csv')
##Figure Comparing Alignment Rate for HISAT2 & kallisto----
data <- read.csv('C:/Users/cnewt/Yang/Dawei_RNA/kallisto_v_hisat.csv')
str(data)
data$Analysis <- as.factor(data$Analysis)
data$Genotype <- as.factor(data$Genotype)
data$Time <- as.factor(data$Time)
data_sum <- ddply(data, c('Analysis', 'Genotype', 'Time'), summarize, Percent = mean(Alignment), STD = sd(Alignment))
str(data_sum)
data_sum$Analysis <- as.factor(data_sum$Analysis)
data_sum$Genotype <- as.factor(data_sum$Genotype)
data_sum$Time <- as.factor(data_sum$Time)
fig <- ggplot(data_sum, aes(x = Time, y = Percent, fill = Analysis)) +
  geom_col(position = 'dodge') +
  facet_wrap(~Genotype) + 
  theme_classic() + 
  ylab('Overall Alignment Rate of Reads (%)') + 
  xlab('Time (hours)') +
  geom_errorbar(aes(ymin=Percent-STD, ymax=Percent+STD, group = Analysis), width = 0.5, position = position_dodge(width = 0.9))
fig
ggsave('KvH.jpeg', device = 'jpeg', scale = 1, width = 6, height = 6, dpi = 'print')
