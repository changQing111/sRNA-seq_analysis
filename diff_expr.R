library(tidyverse)
gene_expr_matrix <- read.csv('./gene_results.csv', header = T)
transcript_expr_matrix <- read.csv('./transcript_results.csv', header = T)
names(gene_expr_matrix) <- c('gene_id','dj_mock_1','dj_mock_2','dj_mock_3','dcl_ab_1', 'dcl_ab_2', 'dcl_ab_3')
rownames(gene_expr_matrix) <- gene_expr_matrix$gene_id
names(transcript_expr_matrix) <- c('transcript_id','dj_mock_1','dj_mock_2','dj_mock_3','dcl_ab_1', 'dcl_ab_2', 'dcl_ab_3')

# DESeq2
library(DESeq2)
# count matrix
cts <- as.matrix(gene_expr_matrix[,2:7]+1) 
cts <- (gene_expr_matrix %>% filter(!(dj_mock_1==0&dj_mock_2==0&dj_mock_3==0&dcl_ab_1==0&dcl_ab_2==0&dcl_ab_3==0)))[,2:7]
# 
condition <- factor(c('dj_mock','dj_mock','dj_mock','dcl_ab', 'dcl_ab', 'dcl_ab'), levels = c('dj_mock', 'dcl_ab'))
# a table of sample information
coldata <- data.frame(row.names = colnames(cts), condition)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="condition_dcl_ab_vs_dj_mock")
res <- res[order(res$padj),]
#diff_gene <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
#diff_expr <- as.data.frame(diff_gene@listData)
#rownames(diff_expr) <- diff_gene@rownames
all_diff_expr <- as.data.frame(res@listData)
rownames(all_diff_expr) <- res@rownames
all_diff_expr <- na.omit(all_diff_expr)
#火山图
all_diff_expr$change <- if (all_diff_expr$padj < 0.05) {
              ifelse (all_diff_expr$log2FoldChange > 1.5,  'Up',
              ifelse (all_diff_expr$log2FoldChange< -1.5,  'Down', 'Stable'))
}
p <- ggplot(na.omit(all_diff_expr[-1,]), aes(x=log2FoldChange, y=-log10(padj))) + geom_point(aes(colour = change), alpha=0.5)
p + geom_vline(xintercept = c(-1.5, 1.5), linetype="dashed", size=0.75) + 
  scale_color_manual(values=c("#00008B","#808080","#DC143C")) +
   theme_bw() + xlim(c(-15, 15)) + ylim(c(0, 20)) + guides(colour=guide_legend(title = NULL))
