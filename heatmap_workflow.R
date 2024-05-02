work <- ("/work/lylab/cjn40747/dawei_RNA/counts")
setwd(work)

library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(tibble)

set.seed(1)

#Data loading and prep----
input <- dir(pattern = '\\.tsv$')
for (file in input) {
  name <- tools::file_path_sans_ext(file)
  assign(name, read.delim(file, header = TRUE, sep = '\t'))
}
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
compiled <- Reduce(function(x,y) merge(x, y, all=TRUE), dfs)
compiled <- column_to_rownames(compiled, 'target_id')

meta <- read.csv('meta.csv')
str(meta)
meta$Genotype <- as.factor(meta$Genotype)
meta$Time <- as.factor(meta$Time)

hour0 <- as.matrix(compiled[, c(1:3, 22:24, 43:45, 64:66)])
min10 <- as.matrix(compiled[, c(4:6, 25:27, 46:48, 67:69)])
hour1 <- as.matrix(compiled[, c(7:9, 28:30, 49:51, 70:72)])
hour24 <- as.matrix(compiled[, c(10:12, 31:33, 52:54, 73:75)])
min30 <- as.matrix(compiled[, c(13:15, 34:36, 55:57, 76:78)])
hour48 <- as.matrix(compiled[, c(16:18, 37:39, 58:60, 79:81)])
hour6 <- as.matrix(compiled[, c(19:21, 40:42, 61:63, 82:84)])

#All treatments----
counts <- DGEList(compiled, samples = meta)
keep_compiled <- filterByExpr(counts) 
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
Group <- factor(paste(meta$Genotype,meta$Time,sep='.'))
meta_compiled <- cbind(meta, Group=Group)
meta_compiled <- model.matrix(~0+Group)
colnames(meta_compiled) <- levels(Group)
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
  levels = meta_compiled)
col10_deg <- glmQLFTest(fit_compiled, contrast = compiled_contrasts[,'col10'])
col10_deg_filter <- topTags(col10_deg)
FDR <- p.adjust(col10_deg$table$PValue, method = 'BH')
sum(FDR < 0.05)
summary(decideTests(col10_deg))
plotMD(col10_deg)