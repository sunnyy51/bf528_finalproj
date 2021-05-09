####Individual Final Project####
#Project 4 Analyst role
library(dplyr)
library(scatterplot3d)
library(ggplot2)
library(Seurat)
library(R.filesets)

#Load data and identify marker genes for each cluster 
ds <- loadRDS("GSM2230760_seurat.rda")
DEGs <- FindAllMarkers(ds, only.pos=TRUE, min.pct = 0.15, logfc.threshold = 0.25)
DEGs %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
head(DEGs, n = 5)

#Filter only significant markers
DEGs_sig<- DEGs[DEGs$p_val_adj<0.05,]

#Save to csv file 
write.csv(DEGs_sig, file="DEGs_sig.csv")

#Clustering 
#find top markers for every cluster 
DEGs <- FindAllMarkers(ds, only.pos=TRUE, min.pct = 0.15, logfc.threshold = 0.25)
top5<-DEGs %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10<-DEGs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- DEGs %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

#Delta/Gamma (SST, PPY)= cluster 0
#Acinar (CPA1) = cluster 3
#Alph (GCG)= cluster 4
#Beta (INS) =  cluster 1, 6
#Cytotoxic T (CD3,CD8) = cluster 9
#Ductal (KRT19) = cluster 5, cluster 10
#Macrophage (CD163,CD68, IgG)= cluster 10, 11, 12
#Vascular (VWF, PECAM1, CD34)= cluster 12
#Epsilon (GHRL) - not found
#Stellate (PDGFRB)- not found
#Mast (TPSAB1, KIT, CPA3)- not found

# Violin Plots using features from Paper 
pdf("violins.pdf")
VlnPlot(ds, features = c("GCG")) # Alpha
VlnPlot(ds, features = c("INS")) # Beta
VlnPlot(ds, features = c("SST")) # Delta
VlnPlot(ds, features = c("KRT19")) # Ductal 
VlnPlot(ds, features = c("PPY")) # Gamma
VlnPlot(ds, features = c("CPA1")) #Acinar
VlnPlot(ds, features = c("PDGFRB")) # Stellate
VlnPlot(ds, features = c("VWF", "PECAM1","CD34")) # Vascular
VlnPlot(ds, features = c("CD68", "CD163")) # Macrophage
dev.off()

#Assign cell type identity to clusters and plot UMAP
new.cluster.ids <- c("Alpha", "Beta", "Delta", "Gamma", "Macrophage", "Ductal", "Acinal", "Vascular")
names(new.cluster.ids) <- levels(ds)
ds <- RenameIdents(ds, new.cluster.ids)
DimPlot(ds, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
png("UMAP_plot.png")

#Heatmap to visualize top marker genes per cluster 
top10 <- DEGs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(ds, features = top10$gene) + NoLegend()

#Identify novel marker genes
novel_mg <- FindAllMarkers(ds, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.50,verbose = TRUE,pseudocount.use = 1)
novel_list <- novel_marker_genes %>% group_by(cluster) %>% top_n(n=1, wt=avg_log2FC)
write.csv(novel_list,"novel_marker.csv")