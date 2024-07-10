work <- ("/work/lylab/cjn40747/dawei_RNA/analysis4")
setwd(work)
library(edgeR)
library(dplyr)
library(tximport)
library(ensembldb)
library(GenomicFeatures)

set.seed(1)

samples <- read.csv('h5_files_meta.csv')
rownames(samples) <- samples$sample
files <- file.path('/work/lylab/cjn40747/dawei_RNA/h5_files', samples$sample, 'abundance.h5')
txdb <- makeTxDbFromGFF('/work/lylab/cjn40747/dawei_RNA/TAIR10_DNA.gtf')
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

#10 min----
A4_col10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col10'])
A4_col10_deg_list <- topTags(A4_col10_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_col10_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')
A4_col10_deg_df <- A4_col10_deg_list$table
A4_10_deg_col_df <- A4_col10_deg_df[A4_col10_deg_df$logFC > 1 | A4_col10_deg_df$logFC < -1, ]

A4_sid10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid10'])
A4_sid10_deg_list <- topTags(A4_sid10_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_sid10_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
A4_sid10_deg_df <- A4_sid10_deg_list$table
A4_10_deg_sid_df <- A4_sid10_deg_df[A4_sid10_deg_df$logFC > 1 | A4_sid10_deg_df$logFC < -1, ]

A4_tga10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga10'])
A4_tga10_deg_list <- topTags(A4_tga10_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A1_tga10_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')
A4_tga10_deg_df <- A4_tga10_deg_list$table
A4_10_deg_tga_df <- A4_tga10_deg_df[A4_tga10_deg_df$logFC > 1 | A4_tga10_deg_df$logFC < -1, ]

A4_10_deg_colsid_df <- merge(A4_10_deg_col_df, A4_10_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_10_deg_colsid_notga_df <- anti_join(A4_10_deg_colsid_df, A4_10_deg_tga_df, by = 'Gene')
names(A4_10_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_10_deg_colsid_notga_df, 'A4_10_deg_df.csv')

#30 min----
A4_col30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col30'])
A4_col30_deg_list <- topTags(A4_col30_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_col30_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')
A4_col30_deg_df <- A4_col30_deg_list$table
A4_30_deg_col_df <- A4_col30_deg_df[A4_col30_deg_df$logFC > 1 | A4_col30_deg_df$logFC < -1, ]

A4_sid30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid30'])
A4_sid30_deg_list <- topTags(A4_sid30_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_sid30_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
A4_sid30_deg_df <- A4_sid30_deg_list$table
A4_30_deg_sid_df <- A4_sid30_deg_df[A4_sid30_deg_df$logFC > 1 | A4_sid30_deg_df$logFC < -1, ]

A4_tga30_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga30'])
A4_tga30_deg_list <- topTags(A4_tga30_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A1_tga30_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')
A4_tga30_deg_df <- A4_tga30_deg_list$table
A4_30_deg_tga_df <- A4_tga30_deg_df[A4_tga30_deg_df$logFC > 1 | A4_tga30_deg_df$logFC < -1, ]

A4_30_deg_colsid_df <- merge(A4_30_deg_col_df, A4_30_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_30_deg_colsid_notga_df <- anti_join(A4_30_deg_colsid_df, A4_30_deg_tga_df, by = 'Gene')
names(A4_30_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_30_deg_colsid_notga_df, 'A4_30_deg_df.csv')

#1 hour----
A4_col1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col1'])
A4_col1_deg_list <- topTags(A4_col1_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_col1_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')
A4_col1_deg_df <- A4_col1_deg_list$table
A4_1_deg_col_df <- A4_col1_deg_df[A4_col1_deg_df$logFC > 1 | A4_col1_deg_df$logFC < -1, ]

A4_sid1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid1'])
A4_sid1_deg_list <- topTags(A4_sid1_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_sid1_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
A4_sid1_deg_df <- A4_sid1_deg_list$table
A4_1_deg_sid_df <- A4_sid1_deg_df[A4_sid1_deg_df$logFC > 1 | A4_sid1_deg_df$logFC < -1, ]

A4_tga1_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga1'])
A4_tga1_deg_list <- topTags(A4_tga1_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A1_tga1_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')
A4_tga1_deg_df <- A4_tga1_deg_list$table
A4_1_deg_tga_df <- A4_tga1_deg_df[A4_tga1_deg_df$logFC > 1 | A4_tga1_deg_df$logFC < -1, ]

A4_1_deg_colsid_df <- merge(A4_1_deg_col_df, A4_1_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_1_deg_colsid_notga_df <- anti_join(A4_1_deg_colsid_df, A4_1_deg_tga_df, by = 'Gene')
names(A4_1_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_1_deg_colsid_notga_df, 'A4_1_deg_df.csv')

#6 hour----
A4_col6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col6'])
A4_col6_deg_list <- topTags(A4_col6_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_col6_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')
A4_col6_deg_df <- A4_col6_deg_list$table
A4_6_deg_col_df <- A4_col6_deg_df[A4_col6_deg_df$logFC > 1 | A4_col6_deg_df$logFC < -1, ]

A4_sid6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid6'])
A4_sid6_deg_list <- topTags(A4_sid6_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_sid6_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
A4_sid6_deg_df <- A4_sid6_deg_list$table
A4_6_deg_sid_df <- A4_sid6_deg_df[A4_sid6_deg_df$logFC > 1 | A4_sid6_deg_df$logFC < -1, ]

A4_tga6_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga6'])
A4_tga6_deg_list <- topTags(A4_tga6_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A1_tga6_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')
A4_tga6_deg_df <- A4_tga6_deg_list$table
A4_6_deg_tga_df <- A4_tga6_deg_df[A4_tga6_deg_df$logFC > 1 | A4_tga6_deg_df$logFC < -1, ]

A4_6_deg_colsid_df <- merge(A4_6_deg_col_df, A4_6_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_6_deg_colsid_notga_df <- anti_join(A4_6_deg_colsid_df, A4_6_deg_tga_df, by = 'Gene')
names(A4_6_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_6_deg_colsid_notga_df, 'A4_6_deg_df.csv')

#24 hour----
A4_col24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col24'])
A4_col24_deg_list <- topTags(A4_col24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_col24_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')
A4_col24_deg_df <- A4_col24_deg_list$table
A4_24_deg_col_df <- A4_col24_deg_df[A4_col24_deg_df$logFC > 1 | A4_col24_deg_df$logFC < -1, ]

A4_sid24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid24'])
A4_sid24_deg_list <- topTags(A4_sid24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_sid24_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
A4_sid24_deg_df <- A4_sid24_deg_list$table
A4_24_deg_sid_df <- A4_sid24_deg_df[A4_sid24_deg_df$logFC > 1 | A4_sid24_deg_df$logFC < -1, ]

A4_tga24_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga24'])
A4_tga24_deg_list <- topTags(A4_tga24_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A1_tga24_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')
A4_tga24_deg_df <- A4_tga24_deg_list$table
A4_24_deg_tga_df <- A4_tga24_deg_df[A4_tga24_deg_df$logFC > 1 | A4_tga24_deg_df$logFC < -1, ]

A4_24_deg_colsid_df <- merge(A4_24_deg_col_df, A4_24_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_24_deg_colsid_notga_df <- anti_join(A4_24_deg_colsid_df, A4_24_deg_tga_df, by = 'Gene')
names(A4_24_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_24_deg_colsid_notga_df, 'A4_24_deg_df.csv')

#48 hour----
A4_col48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col48'])
A4_col48_deg_list <- topTags(A4_col48_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_col48_deg_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0')
A4_col48_deg_df <- A4_col48_deg_list$table
A4_48_deg_col_df <- A4_col48_deg_df[A4_col48_deg_df$logFC > 1 | A4_col48_deg_df$logFC < -1, ]

A4_sid48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'sid48'])
A4_sid48_deg_list <- topTags(A4_sid48_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A4_sid48_deg_df) <- c('Gene', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
A4_sid48_deg_df <- A4_sid48_deg_list$table
A4_48_deg_sid_df <- A4_sid48_deg_df[A4_sid48_deg_df$logFC > 1 | A4_sid48_deg_df$logFC < -1, ]

A4_tga48_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'tga48'])
A4_tga48_deg_list <- topTags(A4_tga48_deg, sort.by = 'logFC', n = 'Inf', p.value = 0.05)
names(A1_tga48_deg_df) <- c('Gene', 'logFC_tga256', 'logCPM_tga256', 'F_tga256', 'Pvalue_tga256', 'FDR_tga256')
A4_tga48_deg_df <- A4_tga48_deg_list$table
A4_48_deg_tga_df <- A4_tga48_deg_df[A4_tga48_deg_df$logFC > 1 | A4_tga48_deg_df$logFC < -1, ]

A4_48_deg_colsid_df <- merge(A4_48_deg_col_df, A4_48_deg_sid_df, by.x = 'Gene', by.y = 'Gene')
A4_48_deg_colsid_notga_df <- anti_join(A4_48_deg_colsid_df, A4_48_deg_tga_df, by = 'Gene')
names(A4_48_deg_colsid_notga_df) <- c('Gene', 'logFC_Col0', 'logCPM_Col_0', 'F_Col0', 'Pvalue_Col0', 'FDR_Col0', 'logFC_sid2', 'logCPM_sid2', 'F_sid2', 'Pvalue_sid2', 'FDR_sid2')
write.csv(A4_48_deg_colsid_notga_df, 'A4_48_deg_df.csv')