library(clusterProfiler)

# data process
df <- read.table("Intersect_gene_ENSG_of_SV.bed", header = F)
colnames(df) <- "ENSEMBL"
gene.df = bitr(df$ENSEMBL, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

# plot
pdf(file = "GO_enrichment_result.pdf", width=5,height=8)
dotplot(ego, showCategory=20)
dev.off()