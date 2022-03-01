# Process sc data, UMAP clustering, DEG analysis
# Basic seurat libs
library(dplyr)
library(Seurat)
library(patchwork)
# plot libs
library(tidyverse)
library(Matrix)
library(cowplot)
library(viridis)
theme_set(theme_cowplot())

dir.create("figures")
Project.name <- "s4T48"
# Load the dataset
data <- Read10X(data.dir = "/zfs/lasernode/feltuslab/gaoyy/mtr_sc_RNASeq/Time_48h/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
seurat_obj <- CreateSeuratObject(
  counts = data,
  project = Project.name,
  min.cells =3,
  min.features=200
)

# plot distributions of QC metrics
png(paste0('figures/',Project.name,'_basic_qc.png'), width=10, height=5, res=200, units='in')
vln <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2)#, pt.size=0)
print(vln)
dev.off()

# apply filter
# nCount_RNA > 200
# nCount_RNA < 2500
seurat_obj <- subset(seurat_obj, nCount_RNA > 200 & nCount_RNA < 2500 )

# log normalize data
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# scale data:
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

#Identification of highly variable features (Feature Selection)
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

png(paste0('figures/',Project.name,'_variable_features.png'), width=5, height=5, res=200, units='in')
p <- LabelPoints(
  VariableFeaturePlot(seurat_obj),
  points = head(VariableFeatures(seurat_obj),10),
  repel = TRUE
) + theme(legend.position="bottom")
print(p)
dev.off()

# Perform linear dimensional reduction
seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(object = seurat_obj)
)
# Print out top 5 PCA
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize the TOP3 PCA
p <- VizDimLoadings(seurat_obj, dims = 1:3, reduction = "pca", ncol=3)

png(paste0('figures/',Project.name,'_pca_loadings.png'), width=15, height=5, res=200, units='in')
print(p)
dev.off()

# PCA scatter plot 
p <- DimPlot(seurat_obj, reduction = "pca")

png(paste0('figures/',Project.name,'_pca_scatter.png'), width=7, height=6, res=200, units='in')
print(p)
dev.off()

# plot heatmaps for first 16 PCs
png(paste0('figures/',Project.name,'_pca_heatmap.png'), width=10, height=10, res=200, units='in')
p <- DimHeatmap(seurat_obj, dims = 1:16, cells = 500, balanced = TRUE, ncol=4)
print(p)
dev.off()

# Plot elbow plot
png(paste0('figures/',Project.name,'_pca_elbow.png'), width=6, height=3, res=200, units='in')
p <- ElbowPlot(seurat_obj, ndims = 50)
print(p)
dev.off()

# Cluster the cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.125)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot UMAP Cluster
umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()
)

png(paste0('figures/',Project.name,'_inital_umap_clusters.png'), width=7, height=7, res=200, units='in')
p <- DimPlot(seurat_obj, reduction = "umap", label=TRUE) +
  umap_theme + NoLegend() + ggtitle('UMAP colored by seurat clusters')
print(p)
dev.off()

# Remove Cluster 0 from the UMAP plot
initial.cluster <- levels(x = seurat_obj)
indices <- c(1)
new.cluster <- initial.cluster[-indices]

# Run UMAP again after excluding Cluster 0
png(paste0('figures/',Project.name,'_sorted_umap_clusters.png'), width=7, height=7, res=200, units='in')
seurat_obj <- subset(seurat_obj,idents = new.cluster)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
p <- DimPlot(seurat_obj, reduction = "umap", label=TRUE) +
  umap_theme + NoLegend() + ggtitle('UMAP colored by seurat clusters')
print(p)
dev.off()

#Finding differentially expressed features (cluster biomarkers)
seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Find the Top3 Markers for each cluster
Top3_markers_df <- seurat_obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

Cluster_markers <- Top3_markers_df$gene

plot_list <- FeaturePlot(seurat_obj, features = Cluster_markers,combine=FALSE, cols=viridis(256))

# apply theme to each feature plot
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + umap_theme + NoLegend()
}

png(paste0('figures/basic_',Project.name,'_Top3marker_featurePlot.png'), width=16, height=10, units='in', res=200)
p<- CombinePlots(plot_list)
print(p)
dev.off()
# Save the marker list 
write.csv(seurat_obj.markers, file=paste0(Project.name,'_cluster_markers.csv'), row.names=T)


# plot the number of DEGs per cluster:
df <- as.data.frame(rev(table(seurat_obj.markers$cluster)))
colnames(df) <- c('cluster', 'n_DEGs')
# bar plot of the number of DEGs in cluster
p <- ggplot(df, aes(y=n_DEGs, x=reorder(cluster, -n_DEGs), fill=cluster)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  NoLegend() + RotatedAxis() +
  ylab(expression(italic(N)[DEGs])) + xlab('') +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="lightgray", size=0.5),
  )

png(paste0('figures/',Project.name,'_DEGs_barplot.png'), width=9, height=4, res=300, units='in')
print(p)
dev.off()

# plot the top 3 DEGs per cluster as a heatmap:
top_DEGs <- seurat_obj.markers %>%
  group_by(cluster) %>%
  top_n(3, wt=avg_log2FC) %>%
  .$gene

png(paste0('figures/',Project.name,'_DEGs_heatmap.png'), width=10, height=10, res=300, units='in')
p <- DoHeatmap(seurat_obj, features=top_DEGs, group.by='seurat_clusters', label=T) + scale_fill_gradientn(colors=viridis(256)) + NoLegend()
print(p)
dev.off()
