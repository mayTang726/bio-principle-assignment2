# Assignment Script

# Change working directory.
path  = "/Users/stan/Desktop/bio principle/final- report - data analysis" # change this to your own directory
file_name = "brca_tcga_pan_can_atlas_2018.tar.gz"

# extract the files into folders.
untar(file_name)

# change directory to the extracted folders
setwd(paste(getwd() , "/brca_tcga_pan_can_atlas_2018", sep = ""))

# data_clinical_patient.txt, data_mrna_seq_v2_rsem.txt, data_mutations.txt and data_cna.txt
clinical = read.delim("data_clinical_patient.txt")
rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")
# in this assignment we will delete the genes for which there's more than one Hugo Symbol
# These are typically genes with no Hugo Symbol ("" as an entry) or pseudogenes.

# This is more for simplicity.If you keep your analysis would still be correct so no worries.
keep = !duplicated(rnaseq[,1])
rnaseq = rnaseq[keep,]

# set rownames of rnaseq to hugo symbols
rownames(rnaseq)  = rnaseq[,1]

# Read CNA Data
cna = read.delim('data_cna.txt')

# find ERBB2 in cna
erbb2_indx = which(cna[,1] == 'ERBB2')

# Plot histogram to visualize explore the data.
hist(as.numeric(cna[erbb2_indx,-c(1,2)]), main = "Histogram of ERBB2 CNA Levels")
# match patients in rnaseq to patients in cna.
rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), colnames(cna[,-c(1,2)])))

# select only the rna cases which have cna data.
rna_cna_sub = rnaseq[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna
no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(rnaseq[,2+rna_cna_id]), colnames(cna[,-c(1,2)]))) 

# sanity check.This will print an error if the result is not the same.
sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]

# Pre-allocate memory for ERBB2
meta_erbb2 = matrix(0,length(rna_cna_id),1)

for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
}

print(meta_erbb2)

# This are some checks you can do to make sure your code worked.
# There's some more systematic checks you can do. See unit testing.


# simple checks to make sure. 
col_i = colnames(rna_cna_sub)[1]
col_cna = which(colnames(cna)==col_i)

# sanity check
(cna[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]

# see now if a positive meta_erbb2 is amplified.
pos_example = which(meta_erbb2==1)[1]
print(pos_example)
col_i = colnames(rna_cna_sub)[pos_example]
col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]

# botch checks should print true.

# We will add a title to the metadata.
colnames(meta_erbb2) = 'ERBB2Amp'

# transform into integers
rna_cna_sub = round(rna_cna_sub)
print(rna_cna_sub)

# Install DESeq2.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DeSeq2
BiocManager::install("DESeq2")
library(DESeq2)

#normaliza data 
dds <- DESeqDataSetFromMatrix(countData = round(rna_cna_sub),
                              colData = meta_erbb2,
                              design = ~ERBB2Amp)

dds <- DESeq(dds)

# perform differential expression analysis
res <- results(dds)
print(res)

# Summary
summary(res)
rownames(res) = rnaseq[keep,1]

# Significantly Differentially Expressed
## column 解释
# log2FoldChange：对数二倍差异，表示实验组相对于对照组的表达水平的倍数变化。正值表示上调，负值表示下调。
# lfcSE：log2FoldChange 的标准误差，衡量 log2FoldChange 的不确定性
# stat: 统计值，通常是 log2FoldChange 除以其标准误差，用于检验 log2FoldChange 是否显著不同于零
# pvalue: 用于判断差异是否显著。较小的 p-value 表示较显著的差异。
# padj: 校正过的 p-value，通常采用多重检验校正方法（如Benjamini-Hochberg），以控制假阳性率
##
# p-value < 0.05 
signif = which(res$padj<0.05)
print(signif)
deg = res[signif,]
print(deg)

# sort result,  and find top ten differentially expression genes
logFC_initial <- res$log2FoldChange
abs_logFC <- abs(res$log2FoldChange)
print(abs_logFC)
sorted_index <- order(abs_logFC, decreasing = TRUE)
sorted_res <- res[sorted_index, ]
print(sorted_index)
# select top ten sorted data
# 取log2FoldChange的绝对值取top 10 
top_ten <- sorted_res[1:10, ]
print(top_ten)

# Separate deg
# combine with data ---> pvalue < 0.05(significantly differentially expression) to analysis gene expression
dup = deg[deg[,2]>0.,] # up regulated
ddown = deg[deg[,2]<0.,] # down regulated

# plot bar chart show the percentage of dup and ddown in deg
up_ratio <- nrow(dup)/nrow(deg)
down_ratio <- nrow(ddown)/nrow(deg)
label <- c(paste('dup:',round(up_ratio * 100),"%"),paste('ddown:',round(down_ratio * 100),"%")) 
colors <- c("pink", "green")
pie(c(up_ratio, down_ratio),labels = label,col = colors, main = "Proportion of dup and ddown in deg")

# For Pathway Enrichment we need Entrez IDs
entrez_ids = rnaseq[keep,2]
entrez_all = entrez_ids[signif]
print(entrez_all)
entrez_up = entrez_all[signif[deg[,2]>0.]]
entrez_down = entrez_all[signif[deg[,2]<0.]]
print(entrez_up)
print(entrez_down)

# Pathway Enrichment
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
library(clusterProfiler)
library(pathview)
# Do a KEGG pathway over-representation analysis
all_paths = enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.05)

# plot different gene expression in hippo signaling pathway 
data <- c(logFC_initial)
df <- data.frame(entrezID = data)
pathview(gene.data = df, #上面包含行名为entrezID的logFC值的数据框
         pathway.id = "hsa04390", #选择一个KEGG信号通路
         species = "hsa",
         out.suffix = "4")

print(all_paths)
head(all_paths,10)
dotplot(all_paths)
barplot(all_paths)

# Optionally you can divide between up and down.
# Both options are Ok for the assignment.
up_paths = enrichKEGG(gene = entrez_up, organism = 'hsa', pvalueCutoff = 0.05)
print(up_paths)
# dotplot(up_paths)
head(up_paths) #0 row ????
# dotplot(up_paths)
# barplot(up_paths)

down_paths = enrichKEGG(gene = entrez_down, organism = 'hsa', pvalueCutoff = 0.05)
print(down_paths)
head(down_paths) # 0 row ???
# dotplot(down_paths)
# barplot(down_paths)

# Transform the data to visualize
rld <- vst(dds, blind=FALSE)

# Do Principal Components Analysis  pca 分析
pc = prcomp(assay(rld))

# extract first principal component 
loadings_pc1 <- pc$rotation[, 2]
# combine the results with initial variation
loadings_data <- data.frame(Variable = colnames(assay(rld)), Loadings_PC1 = loadings_pc1)
# sort by absolute value
loadings_data <- loadings_data[order(abs(loadings_data$Loadings_PC1), decreasing = TRUE), ]
print(head(loadings_data))

# Plot 
# pc1 means amplified 
# red dot identify amplified genes, black identify not amplofied genes
plot(pc$rotation[,1], pc$rotation[,2], col = 1+(meta_erbb2), pch = 19)

