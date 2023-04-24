library(DropletUtils)
library(Seurat)
library(ggplot2)
library(scattermore)
library(patchwork)
library(ggridges)
library(ggrepel)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette = rev(cbPalette)
cbPalette = c("#000000", cbPalette)


## ==== Import datasets

importData = function(datadir){
  
  temp = list()
  
  # Import
  sce_filtered = read10xCounts(paste0(datadir,'/filtered_feature_bc_matrix'))
  colnames(sce_filtered) = colData(sce_filtered)$Barcode
  
  sce_raw = read10xCounts(paste0(datadir,'/raw_feature_bc_matrix'))
  colnames(sce_raw) = colData(sce_raw)$Barcode
  
  
  gex_filtered = sce_filtered[rowData(sce_filtered)$Type == 'Gene Expression',]
  gex_raw = sce_raw[rowData(sce_raw)$Type == 'Gene Expression',]
  
  adt_filtered = sce_filtered[rowData(sce_filtered)$Type == 'Antibody Capture',]
  adt_raw = sce_raw[rowData(sce_raw)$Type == 'Antibody Capture',]
  
  
  
  # Filter out all 0 droplets
  df = data.frame(rna_size = colSums2(counts(gex_raw)),
                  prot_size = colSums2(counts(adt_raw)))
  
  rownames(df) = colnames(gex_raw)
  
  df$emptyDrop_rna = ifelse(colnames(gex_raw) %in% colnames(sce_filtered), "cell", "background")
  df$emptyDrop_rna = as.factor(df$emptyDrop_rna)
  
  non_zero_drop = df[df$rna_size > 0 & df$prot_size > 0,]
  non_zero_drop$rna_size_log = log1p(non_zero_drop$rna_size)
  non_zero_drop$prot_size_log = log1p(non_zero_drop$prot_size)
  
  
  temp[['adt_filtered']] = adt_filtered
  temp[['adt_raw']] = adt_raw
  temp[['gex_filtered']] = gex_filtered
  temp[['gex_raw']] = gex_raw
  temp[['non_zero_drop']] = non_zero_drop
  
  return(temp)
  
}



datadir_name = 'data_pbmc10khealthdonor'
datadir_10k = paste0('/rprojectnb2/camplab/home/yin/poisson/', datadir_name)
pbmc10k = importData(datadir_10k)

datadir_name = 'data_pbmc5knextgem'
datadir_5k = paste0('/rprojectnb2/camplab/home/yin/poisson/', datadir_name)
pbmc5k = importData(datadir_5k)

datadir_name = 'data_malt10k'
datadir_malt = paste0('/rprojectnb2/camplab/home/yin/poisson/', datadir_name)
malt10k = importData(datadir_malt)

zhang = importData('/rprojectnb2/camplab/home/yin/decontX_ADT/zhangetal/young-aged/')




## ==== Clustering

# PBMC10k
non_zero_drop = pbmc10k$non_zero_drop

# Use K-means to cluster empty drop
temp = non_zero_drop[non_zero_drop$emptyDrop_rna == 'background' &
                       (non_zero_drop$rna_size_log > 1 | non_zero_drop$prot_size_log > 1),]

cl = kmeans(temp[,c('rna_size_log', 'prot_size_log')],
            centers = matrix(c(4,0,2,5,3,9), nrow = 3, ncol = 2, byrow = T))

temp$cluster = cl$cluster

# Check
scattermoreplot(temp$rna_size_log,
                temp$prot_size_log,
                col = cbPalette[c(1,2)][as.integer(temp$cluster == 1) + 1],
                cex=0.5,
                xlab = 'RNA Lib Size',
                ylab = 'ADT Lib Size')

non_zero_drop$emptyDrop_high = rownames(non_zero_drop) %in% rownames(temp[temp$cluster==3,])
non_zero_drop$emptyDrop_mid = rownames(non_zero_drop) %in% rownames(temp[temp$cluster==2,])
non_zero_drop$emptyDrop_low = non_zero_drop$emptyDrop_rna == 'background' &
  !non_zero_drop$emptyDrop_high &
  !non_zero_drop$emptyDrop_mid


pbmc10k$non_zero_drop = non_zero_drop




# PBMC5k
non_zero_drop = pbmc5k$non_zero_drop

# Use K-means to cluster empty drop
temp = non_zero_drop[non_zero_drop$emptyDrop_rna == 'background' &
                       (non_zero_drop$rna_size_log > 1 | non_zero_drop$prot_size_log > 1),]

cl = kmeans(temp[,c('rna_size_log', 'prot_size_log')],
            centers = matrix(c(2,2,4,4), nrow = 2, ncol = 2, byrow = T))

temp$cluster = cl$cluster

# Check
scattermoreplot(temp$rna_size_log,
                temp$prot_size_log,
                col = cbPalette[c(1,2)][as.integer(temp$cluster == 2) + 1],
                cex=0.5,
                xlab = 'RNA Lib Size',
                ylab = 'ADT Lib Size')


non_zero_drop$emptyDrop_mid = rownames(non_zero_drop) %in% rownames(temp[temp$cluster==2,])
non_zero_drop$emptyDrop_low = non_zero_drop$emptyDrop_rna == 'background' &
  !non_zero_drop$emptyDrop_mid


pbmc5k$non_zero_drop = non_zero_drop



# MALT10k
non_zero_drop = malt10k$non_zero_drop

# Use K-means to cluster empty drop
temp = non_zero_drop[non_zero_drop$emptyDrop_rna == 'background' &
                       (non_zero_drop$rna_size_log > 1 | non_zero_drop$prot_size_log > 1),]

cl = kmeans(temp[,c('rna_size_log', 'prot_size_log')],
            centers = matrix(c(2,2,3,5,5,7), nrow = 3, ncol = 2, byrow = T))

temp$cluster = cl$cluster

# Check
scattermoreplot(temp$rna_size_log,
                temp$prot_size_log,
                col = cbPalette[c(1,2)][as.integer(temp$cluster == 2) + 1],
                cex=0.5,
                xlab = 'RNA Lib Size',
                ylab = 'ADT Lib Size')

non_zero_drop$emptyDrop_high = rownames(non_zero_drop) %in% rownames(temp[temp$cluster==3,])
non_zero_drop$emptyDrop_mid = rownames(non_zero_drop) %in% rownames(temp[temp$cluster==2,])
non_zero_drop$emptyDrop_low = non_zero_drop$emptyDrop_rna == 'background' &
  !non_zero_drop$emptyDrop_high &
  !non_zero_drop$emptyDrop_mid


malt10k$non_zero_drop = non_zero_drop





# Zhang
non_zero_drop = zhang$non_zero_drop


non_zero_drop$emptyDrop_high = ifelse(non_zero_drop$emptyDrop_rna == 'background' &
                                        non_zero_drop$prot_size_log > 6,
                                      T, F)
non_zero_drop$emptyDrop_mid = ifelse(non_zero_drop$emptyDrop_rna == 'background' &
                                       non_zero_drop$prot_size_log <= 6 &
                                       non_zero_drop$prot_size_log > 3.2,
                                     T, F)
non_zero_drop$emptyDrop_low = ifelse(non_zero_drop$emptyDrop_rna == 'background' &
                                       non_zero_drop$prot_size_log <= 3.2,
                                     T, F)

zhang$non_zero_drop = non_zero_drop







## ====================
## ====== Figure
## ====================



## ==== Correlation plot

calc_profile = function(sce){
  profile = rowSums2(counts(sce))
  profile = profile/sum(profile)
}


dat4plotting = function(adt_raw, non_zero_drop){
  
  low_bg = adt_raw[, colnames(adt_raw) %in% rownames(non_zero_drop[non_zero_drop$emptyDrop_low,])]
  cell = adt_raw[, colnames(adt_raw) %in% rownames(non_zero_drop[non_zero_drop$emptyDrop_rna=='cell',])]
  
  
  low_bg_profile = calc_profile(low_bg)
  cell_profile = calc_profile(cell)
  
  df = data.frame(low = low_bg_profile,
                  cell = cell_profile,
                  row.names = rownames(adt_raw))
  
}


# PBMC10k
adt_raw = pbmc10k$adt_raw
non_zero_drop = pbmc10k$non_zero_drop

df = dat4plotting(adt_raw, non_zero_drop)
write.csv(df, 'ambient_data_pbmc10khealthdonor.csv')

S1A = ggplot(df, aes(low, cell)) +
  geom_point(shape = 1, size = 2) +
  annotate("text", x=0.25, y=0.01, label=paste0("R = ", round(cor(df$low, df$cell),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  geom_text_repel(label = rownames(df)) +
  ggtitle(paste0('PBMC10k: Cell Average Profile vs. Ambient Profile (ADT)')) + 
  labs(x = 'Ambient Profile', y = 'Cell Average Profile', tag = "A") +
  coord_cartesian(xlim = c(0,0.3), ylim = c(0, 0.3)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold",
                                size = 15))


gex_raw = pbmc10k$gex_raw
non_zero_drop = pbmc10k$non_zero_drop

df = dat4plotting(gex_raw, non_zero_drop)

S1B = ggplot(df, aes(low, cell)) +
  geom_scattermore(pointsize = 2,
                   pixels=c(1000,1000)) +
  annotate("text", x=0.025, y=0.001, label=paste0("R = ", round(cor(df$low, df$cell),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  # geom_text_repel(label = rownames(df)) +
  ggtitle(paste0('PBMC10k: Cell Average Profile vs. Ambient Profile (GEX)')) + 
  labs(x = 'Ambient Profile', y = 'Cell Average Profile', tag = "B") +
  coord_cartesian(xlim = c(0,0.03), ylim = c(0, 0.03)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold",
                                size = 15))


# PBMC5k
adt_raw = pbmc5k$adt_raw
non_zero_drop = pbmc5k$non_zero_drop

df = dat4plotting(adt_raw, non_zero_drop)
write.csv(df, 'ambient_data_pbmc5knextgem.csv')

S1C = ggplot(df, aes(low, cell)) +
  geom_point(shape = 1, size = 2) +
  annotate("text", x=0.25, y=0.01, label=paste0("R = ", round(cor(df$low, df$cell),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  geom_text_repel(label = rownames(df)) +
  ggtitle(paste0('PBMC5k: Cell Average Profile vs. Ambient Profile (ADT)')) + 
  labs(x = 'Ambient Profile', y = 'Cell Average Profile', tag = "C") +
  coord_cartesian(xlim = c(0,0.3), ylim = c(0, 0.3)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold",
                                size = 15))


gex_raw = pbmc5k$gex_raw
non_zero_drop = pbmc5k$non_zero_drop

df = dat4plotting(gex_raw, non_zero_drop)

S1D = ggplot(df, aes(low, cell)) +
  geom_scattermore(pointsize = 2,
                   pixels=c(1000,1000)) +
  annotate("text", x=0.025, y=0.001, label=paste0("R = ", round(cor(df$low, df$cell),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  # geom_text_repel(label = rownames(df)) +
  ggtitle(paste0('PBMC5k: Cell Average Profile vs. Ambient Profile (GEX)')) + 
  labs(x = 'Ambient Profile', y = 'Cell Average Profile', tag = "D") +
  coord_cartesian(xlim = c(0,0.03), ylim = c(0, 0.03)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold",
                                size = 15))



# MALT10k
adt_raw = malt10k$adt_raw
non_zero_drop = malt10k$non_zero_drop

df = dat4plotting(adt_raw, non_zero_drop)
write.csv(df, 'ambient_data_malt10k.csv')

S1E = ggplot(df, aes(low, cell)) +
  geom_point(shape = 1, size = 2) +
  annotate("text", x=0.25, y=0.01, label=paste0("R = ", round(cor(df$low, df$cell),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  geom_text_repel(label = rownames(df)) +
  ggtitle(paste0('MALT10k: Cell Average Profile vs. Ambient Profile (ADT)')) + 
  labs(x = 'Ambient Profile', y = 'Cell Average Profile', tag = "E") +
  coord_cartesian(xlim = c(0,0.3), ylim = c(0, 0.3)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold",
                                size = 15))


gex_raw = malt10k$gex_raw
non_zero_drop = malt10k$non_zero_drop

df = dat4plotting(gex_raw, non_zero_drop)

S1F = ggplot(df, aes(low, cell)) +
  geom_scattermore(pointsize = 2,
                   pixels=c(1000,1000)) +
  annotate("text", x=0.025, y=0.001, label=paste0("R = ", round(cor(df$low, df$cell),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  # geom_text_repel(label = rownames(df)) +
  ggtitle(paste0('MALT10k: Cell Average Profile vs. Ambient Profile (GEX)')) + 
  labs(x = 'Ambient Profile', y = 'Cell Average Profile', tag = "F") +
  coord_cartesian(xlim = c(0,0.03), ylim = c(0, 0.03)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold",
                                size = 15))



# Zhang
adt_raw = zhang$adt_raw
non_zero_drop = zhang$non_zero_drop

df = dat4plotting(adt_raw, non_zero_drop)

S1G = ggplot(df, aes(low, cell)) +
  geom_point(shape = 1, size = 2) +
  annotate("text", x=0.25, y=0.01, label=paste0("R = ", round(cor(df$low, df$cell),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  geom_text_repel(label = rownames(df)) +
  ggtitle(paste0('Golomb: Cell Average Profile vs. Ambient Profile (ADT)')) + 
  labs(x = 'Ambient Profile', y = 'Cell Average Profile', tag = "G") +
  coord_cartesian(xlim = c(0,0.3), ylim = c(0, 0.3)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold",
                                size = 15))


gex_raw = zhang$gex_raw
non_zero_drop = zhang$non_zero_drop

df = dat4plotting(gex_raw, non_zero_drop)

S1H = ggplot(df, aes(low, cell)) +
  geom_scattermore(pointsize = 2,
                   pixels=c(1000,1000)) +
  annotate("text", x=0.025, y=0.001, label=paste0("R = ", round(cor(df$low, df$cell),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  # geom_text_repel(label = rownames(df)) +
  ggtitle(paste0('Golomb: Cell Average Profile vs. Ambient Profile (GEX)')) + 
  labs(x = 'Ambient Profile', y = 'Cell Average Profile', tag = "H") +
  coord_cartesian(xlim = c(0,0.03), ylim = c(0, 0.03)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold",
                                size = 15))




l = list()
l[[1]] = S1A
l[[2]] = S1B
l[[3]] = S1C
l[[4]] = S1D
l[[5]] = S1E
l[[6]] = S1F
l[[7]] = S1G
l[[8]] = S1H


design = "
AABB
CCDD
EEFF
GGHH
"

S1 = wrap_plots(l) + 
  plot_layout(design = design)

ggsave("FigS1.pdf",
       S1,
       width = 8.5*1.3,
       height = 11*1.3,
       units = "in")







## ==== violin plot of zhang dataset

adt_raw = zhang$adt_raw
non_zero_drop = zhang$non_zero_drop

low_bg = adt_raw[, colnames(adt_raw) %in% rownames(non_zero_drop[non_zero_drop$emptyDrop_low,])]
cell = adt_raw[, colnames(adt_raw) %in% rownames(non_zero_drop[non_zero_drop$emptyDrop_rna=='cell',])]


df_cell = reshape2::melt(as.matrix(t(counts(cell))))

S12A = ggplot(df_cell, aes(Var2, value, fill = Var2)) +
  geom_violin(scale = "width") +
  coord_cartesian(ylim = c(0,500)) +
  labs(x="", tag = "A")+
  ggtitle("Golomb: ADT Level in Cell Droplets")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        plot.tag = element_text(face = "bold",
                                size = 15),
        legend.position = "none")




## UMAP and feature plot for Zhang dataset

adt_raw = zhang$adt_raw
non_zero_drop = zhang$non_zero_drop

cells = non_zero_drop[non_zero_drop$emptyDrop_rna == 'cell',]
cells = rownames(cells)

# Seurat to get clusters
temp = adt_raw[,colnames(adt_raw) %in% cells]
adt_seurat = Seurat::CreateSeuratObject(counts(temp),
                                        assay = 'ADT')
adt_seurat <- Seurat::NormalizeData(adt_seurat, 
                                    normalization.method = "CLR",
                                    margin = 2)
adt_seurat <- ScaleData(adt_seurat, 
                        assay = "ADT", 
                        verbose = FALSE)
adt_seurat <- RunPCA(adt_seurat,
                     assay = "ADT",
                     features = rownames(adt_seurat),
                     reduction.name = "pca_adt",
                     reduction.key = "pcaadt_", 
                     verbose = FALSE,
                     npcs = 20)


adt_seurat <- FindNeighbors(adt_seurat,
                            dims = 1:10,
                            assay = "ADT",
                            reduction = "pca_adt",
                            verbose = FALSE)
adt_seurat <- FindClusters(adt_seurat,
                           resolution = 0.2,
                           verbose = FALSE)

adt_seurat <- RunUMAP(adt_seurat,
                      dims = 1:10,
                      assay = "ADT",
                      reduction = "pca_adt",
                      reduction.name = "adtUMAP",
                      verbose = FALSE)

adt_seurat[["adtClusterID"]] <- Idents(adt_seurat)


S12B = DimPlot(adt_seurat, 
              reduction = "adtUMAP",
              label = TRUE,
              repel = TRUE,
              group.by = "adtClusterID") +
  NoLegend() + 
  labs(tag = "B") +
  theme(plot.title = element_blank(),
        plot.tag = element_text())



features = c("CITE-CD11b", "CITE-CD45", "CITE-CX3CR1", "CITE-Ly6C", "CITE-I-A-I-E", "CITE-CD24")
S12C = FeaturePlot(adt_seurat, features = features, combine = FALSE)

for(i in 1:length(S12C)) {
  S12C[[i]] <- S12C[[i]] + NoLegend() + NoAxes()
}

S12C = cowplot::plot_grid(plotlist = S12C) +
  labs(tag = "C")


l = list()
l[[1]] = S12A
l[[2]] = S12B
l[[3]] = S12C


design = "
AAA
BCC
BCC
"

S12 = wrap_plots(l) + 
  plot_layout(design = design)

ggsave("FigS12.pdf",
       S12,
       width = 8.5,
       height = 11,
       units = "in")


