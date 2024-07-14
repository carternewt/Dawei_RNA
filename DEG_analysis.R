setwd("C:/Users/cnewt/Yang/Dawei_RNA/counts") #For local computer

library(edgeR)
library(tibble)
library(dplyr)
library(tximport)
library(rhdf5)
library(ensembldb)
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
counts <- DGEList(txi.kallisto$counts, samples = meta)
keep_compiled <- filterByExpr(counts, meta_compiled) 
counts <- counts[keep_compiled, , keep.lib.sizes=FALSE]
counts <- calcNormFactors(counts, method='TMM')
counts <- normLibSizes(counts)
counts_matrix <- as.matrix(counts$counts)

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