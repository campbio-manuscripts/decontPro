library(DropletUtils)
library(ggplot2)
library(scattermore)
library(patchwork)
library(Seurat)
library(ggridges)
library(ggrepel)
library(cowplot)
library(ClusterR)
library(dplyr)


set.seed(123)


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette = rev(cbPalette)
cbPalette = c("#000000", cbPalette)



## ======== Empty drop analysis

datadir_name = 'data_pbmc10khealthdonor'
datadir = paste0('/rprojectnb2/camplab/home/yin/poisson/', datadir_name)

datadir_name = 'data_pbmc5knextgem'
datadir = paste0('/rprojectnb2/camplab/home/yin/poisson/', datadir_name)

datadir_name = 'data_malt10k'
datadir = paste0('/rprojectnb2/camplab/home/yin/poisson/', datadir_name)

datadir_name = 'data_Golomb'
datadir = '/rprojectnb2/camplab/home/yin/decontX_ADT/zhangetal/young-aged/'






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




# Check by plotting
scattermoreplot(log1p(non_zero_drop$rna_size),
                log1p(non_zero_drop$prot_size),
                col = as.integer(non_zero_drop$emptyDrop_rna),
                cex=0.5,
                xlab = 'RNA Lib Size',
                ylab = 'ADT Lib Size')



non_zero_drop$rna_size_log = log1p(non_zero_drop$rna_size)
non_zero_drop$prot_size_log = log1p(non_zero_drop$prot_size)




# For PBMC10k

# Use K-means to cluster empty drop
temp = non_zero_drop[non_zero_drop$emptyDrop_rna == 'background' &
                     (non_zero_drop$rna_size_log > 1 | non_zero_drop$prot_size_log > 1),]

cl = kmeans(temp[,c('rna_size_log', 'prot_size_log')],
            centers = matrix(c(4,0,2,5,3,9), nrow = 3, ncol = 2, byrow = T))

temp$cluster = cl$cluster

# Check
scattermoreplot(temp$rna_size_log,
                temp$prot_size_log,
                col = cbPalette[c(1,2)][as.integer(temp$cluster == 3) + 1],
                cex=0.5,
                xlab = 'RNA Lib Size',
                ylab = 'ADT Lib Size')



non_zero_drop$emptyDrop_high = rownames(non_zero_drop) %in% rownames(temp[temp$cluster==3,])
non_zero_drop$emptyDrop_mid = rownames(non_zero_drop) %in% rownames(temp[temp$cluster==2,])
non_zero_drop$emptyDrop_low = non_zero_drop$emptyDrop_rna == 'background' &
                              !non_zero_drop$emptyDrop_high &
                              !non_zero_drop$emptyDrop_mid




# # For PBMC5k

# Use K-means to cluster empty drop
temp = non_zero_drop[non_zero_drop$emptyDrop_rna == 'background' &
                       (non_zero_drop$rna_size_log > 1 | non_zero_drop$prot_size_log > 1),]

cl = kmeans(temp[,c('rna_size_log', 'prot_size_log')],
            centers = matrix(c(2,2,4,4), nrow = 2, ncol = 2, byrow = T))

temp$cluster = cl$cluster

# Check
scattermoreplot(temp$rna_size_log,
                temp$prot_size_log,
                col = cbPalette[c(1,2)][as.integer(temp$cluster == 3) + 1],
                cex=0.5,
                xlab = 'RNA Lib Size',
                ylab = 'ADT Lib Size')



non_zero_drop$emptyDrop_mid = rownames(non_zero_drop) %in% rownames(temp[temp$cluster==2,])
non_zero_drop$emptyDrop_low = non_zero_drop$emptyDrop_rna == 'background' &
  !non_zero_drop$emptyDrop_mid






# # For MALT10k

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




# For Golomb

# Cluster in empty drop
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







# Label clusters by numbers
all_drop_label = rep(0, nrow(non_zero_drop))
all_drop_label[non_zero_drop$emptyDrop_rna == 'background'] = 1
all_drop_label[non_zero_drop$emptyDrop_rna == 'cell'] = 2
all_drop_label[non_zero_drop$emptyDrop_high] = 3
all_drop_label[non_zero_drop$emptyDrop_mid] = 4
all_drop_label[non_zero_drop$emptyDrop_low] = 5
non_zero_drop$all_drop_label = factor(all_drop_label)



F1A = ggplot(non_zero_drop,
             aes(rna_size, prot_size, col = all_drop_label)) +
        geom_scattermore(pointsize = 1.2,
                         pixels=c(1000,1000)) +
        scale_color_manual(labels = c('A: Cell','B: Mislabeled\n    Cell','C: Spongelet','D: Ambient'),
                           values = cbPalette[2:5]) +
        scale_x_continuous(trans = 'log1p', breaks=c(10^(0:5)),
                           labels = function(x) format(x, scientific = FALSE)) +
        scale_y_continuous(trans = 'log1p', breaks=c(10^(0:5)),
                           labels = function(x) format(x, scientific = FALSE)) +
        labs(x = 'RNA Library Size', 
             y = 'ADT Library Size',
             tag = "A") +
        geom_vline(xintercept = 600) +
        annotate("text", x = 100, label = 'Empty\nDroplets', y=90000, size = 3) +
        annotate("text", x = 2200, label = 'Filtered\nDroplets', y=90000, size = 3) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size = 12),
              legend.justification = c("right", "bottom"),
              legend.position = c(0.98, 0.02),
              legend.title = element_blank(),
              legend.background = element_blank(),
              legend.box.background = element_rect(colour = "black"),
              legend.text = element_text(margin = margin(l = -10, unit = "pt")),
              plot.tag = element_text(face = "bold",
                                      size = 15))


## Save F1A of PBMC5K, MALT10K, Golomb for FigS3
save(F1A, file = paste0('FigS3_',datadir_name,'.rdata'))







# ================= Analyzing empty drop cluster

set.seed(123)

cells0 = non_zero_drop[non_zero_drop$emptyDrop_rna == 'cell',]
cells0 = rownames(cells0)

cells1 = non_zero_drop[non_zero_drop$emptyDrop_high,]
cells1_all = rownames(cells1)
cells1_top = rownames(cells1[order(cells1$prot_size, decreasing = T)[1:1000],]) # pick top 1000 drop
cells1 = sample(rownames(cells1), min(1000, nrow(cells1)))


cells2 = non_zero_drop[non_zero_drop$emptyDrop_mid,]
cells2_all = rownames(cells2)
cells2_top = rownames(cells2[order(cells2$prot_size, decreasing = T)[1:1000],]) # pick top 1000 drop
cells2 = sample(rownames(cells2), 1000)

cells3 = non_zero_drop[non_zero_drop$emptyDrop_low,]
cells3_all = rownames(cells3)
cells3_top = rownames(cells3[order(cells3$prot_size, decreasing = T)[1:1000],]) # pick top 1000 drop
cells3 = sample(rownames(cells3), 1000)


cells = c(cells0, cells1, cells2)


# Verify cell selection
par(mfrow = c(1,1))
scattermoreplot(non_zero_drop$rna_size_log,
                non_zero_drop$prot_size_log,
                col = as.integer(rownames(non_zero_drop) %in% cells)+1,
                cex=0.5,
                xlab = 'RNA Lib Size',
                ylab = 'ADT Lib Size')





## ========== Clustering using Seurat

temp = adt_raw[,cells]
adt_seurat = Seurat::CreateSeuratObject(counts(temp),
                                        assay = 'ADT')
adt_seurat <- Seurat::NormalizeData(adt_seurat, 
                                    normalization.method = "CLR",
                                    margin = 2)
adt_seurat <- ScaleData(adt_seurat, assay = "ADT", verbose = FALSE)
adt_seurat <- RunPCA(adt_seurat,
                     assay = "ADT",
                     features = rownames(adt_seurat),
                     reduction.name = "pca_adt",
                     reduction.key = "pcaadt_", 
                     verbose = FALSE,
                     npcs = 10)

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

p1 = DimPlot(adt_seurat, reduction = "adtUMAP", label = TRUE, group.by = "adtClusterID") +
  theme(plot.title = element_blank())

features = c('CD3', # Tcell
             'CD4',
             'CD8a',
             'CD19', # Bcell
             'CD14', # Monocyte
             'CD15',
             'CD16',
             'CD25', 
             'CD56', # NK cell and others
             'CD45RA',
             'CD45RO')


#
p2 = FeaturePlot(adt_seurat, features = features)
p1 + p2 + plot_layout(nrow=1,widths=c(1,3))



# Annotate PBMC10k cell type
cell_type = Idents(adt_seurat)
levels(cell_type) = list(MonoCD14 = '0',
                         CD4TMemory = '1',
                         CD4TNative = '2',
                         ClusterC = '3',
                         CD8TMemory = '4',
                         NK = '5',
                         ClusterB = '6',
                         B = '7',
                         CD8TNaive = '8',
                         Doublets = '9',
                         Doublets2 = '10')
adt_seurat[["cellType"]] = cell_type



# All cell cluster 1 color, high and mid empty their own color

label = rep(2, 11)
label[4] = 4
label[7] = 3

F1B = DimPlot(adt_seurat, 
             reduction = "adtUMAP",
             cols = cbPalette[label],
             label = T,
             # label.size = 3,
             repel = T, 
             group.by = "cellType") +
  labs(tag = "B") +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        legend.position = 'none',
        plot.tag = element_text(face = "bold",
                                size = 15))





## ========== Density Plot

# PBMC10k
proteins = c('CD4', 'CD14', 'CD16', 'CD45RA')


p = list()

for (i in 1:length(proteins)) {
  
  protein = proteins[i]
  
  df_bg1 = data.frame(Spongelet=assay(adt_raw[, colnames(adt_raw) %in% cells2_all],'counts')[protein,])
  df_bg1.m = reshape2::melt(df_bg1)
  
  df_bg2 = data.frame(Ambient=assay(adt_raw[, colnames(adt_raw) %in% cells3_top],'counts')[protein,])
  df_bg2.m = reshape2::melt(df_bg2)
  
  df_decon <- data.frame(Cell=assay(adt_raw[, colnames(adt_raw) %in% cells0],'counts')[protein,])
  df_decon.m = reshape2::melt(df_decon)
  
  
  p1 = ggplot(df_decon.m, aes(value, fill = variable)) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(trans=scales::pseudo_log_trans(), breaks=c(1,5,10^(1:4))) +
    ggtitle(protein) +
    scale_fill_manual(values=cbPalette[2]) +
    theme_classic() +
    theme(axis.text = element_text(angle = 0, size = 11),
          axis.title.x=element_blank())
  
  
  # add background
  ylimit = layer_scales(p1)$y$get_limits()
  ylimit[2] = 1.2
  
  p1 = p1 + 
    geom_line(data = df_bg1.m, aes(value, color = variable, fill=NA), 
              stat = 'density',
              adjust = 2,
              size = 1) + 
    theme(legend.title = element_blank(),
          legend.margin=margin(t = -10))
  
  p1 = p1 + 
    geom_line(data = df_bg2.m, aes(value, color = variable, fill=NA), 
              stat = 'density',
              adjust = 2,
              size = 1) + 
    coord_cartesian(ylim = ylimit) +
    scale_color_manual(values=cbPalette[c(5,4)]) + 
    theme(legend.title = element_blank(),
          legend.margin=margin(t = -10))
  
  
  p[[i]] = p1
}

p[[1]] = p[[1]] + labs(tag = "D") + theme(plot.tag = element_text(face = "bold",
                                                                  size = 15))

F1D = wrap_plots(p, ncol = 2, guides = "collect")+ 
  plot_layout()&theme(legend.position = "bottom", 
                      plot.margin = unit(c(3,3,2,1), "pt"))







## Correlation plot

mid_bg = adt_raw[, colnames(adt_raw) %in% cells2_all]
low_bg = adt_raw[, colnames(adt_raw) %in% cells3_all]
cell = adt_raw[, colnames(adt_raw) %in% c(cells0)]



calc_profile = function(sce){
  profile = rowSums2(counts(sce))
  profile = profile/sum(profile)
}


mid_bg_profile = calc_profile(mid_bg)
low_bg_profile = calc_profile(low_bg)
cell_profile = calc_profile(cell)



# 
df = data.frame(mid = mid_bg_profile,
                low = low_bg_profile,
                cell = cell_profile,
                row.names = rownames(adt_raw))


F1C = ggplot(df, aes(mid, low)) +
  geom_point(shape = 1, size = 2) +
  annotate("text", x=0.25, y=0.01, label=paste0("R = ", round(cor(df$mid, df$low),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  geom_text_repel(label = rownames(df)) +
  ggtitle(paste0('Ambient Profile vs. Spongelet Profile')) + 
  labs(x = 'Spongelet Profile', y = 'Ambient Profile', tags = "C") +
  coord_cartesian(xlim = c(0,0.3), ylim = c(0, 0.3)) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        plot.tag = element_text(face = "bold",
                                size = 15))




# PBMC10k

# read off background limit
bg_limit = list(CD3 = c(5, 200),
                CD4 = c(5, 100),
                CD8a = c(5,300),
                CD14 = c(5,100),
                CD15 = c(30, 300),
                CD16 = c(5,300),
                CD56 = c(5, 100),
                CD19 = c(1, 100),
                CD25 = c(3, 100),
                CD45RA = c(10, 400),
                CD45RO = c(5, 100),
                PD1 = c(5, 100),
                TIGIT = c(1, 20),
                CD127 = c(1, 30),
                IgG2a = c(1,30),
                IgG1 = c(1,30),
                IgG2b = c(1,30))



# MALT10k

# read off background limit
bg_limit = list(CD3 = c(5, 300),
                CD4 = c(5, 300),
                CD8a = c(1,15),
                CD14 = c(1,20),
                CD15 = c(1, 80),
                CD16 = c(1,80),
                CD56 = c(5, 80),
                CD19 = c(5, 100),
                CD25 = c(5, 50),
                CD45RA = c(5, 60),
                CD45RO = c(5, 100),
                PD1 = c(5, 100),
                TIGIT = c(1, 30),
                CD127 = c(1, 30),
                IgG2a = c(1,30),
                IgG1 = c(1,30),
                IgG2b = c(1,30))



# PBMC5k

# read off background limit
bg_limit = list(CD3 = c(1, 80),
                CD4 = c(1, 30),
                CD8a = c(1,30),
                CD11b = c(1, 30),
                CD14 = c(1,60),
                CD15 = c(1, 30),
                CD16 = c(0,10),
                CD19 = c(0, 10),
                CD20 = c(1, 30),
                CD25 = c(1, 30),
                CD27 = c(1, 30),
                CD28 = c(1, 30),
                CD34 = c(0, 5),
                CD45RA = c(1, 70),
                CD45RO = c(1, 20),
                CD56 = c(1, 30),
                CD62L = c(1, 15),
                CD69 = c(1, 80),
                CD80 = c(0, 10),
                CD86 = c(0, 20),
                CD127 = c(1, 15),
                CD137 = c(1, 10),
                CD197 = c(1, 20),
                CD274 = c(1, 10),
                CD278 = c(1, 15),
                CD335 = c(1, 10),
                PD1 = c(1, 10),
                HLADR = c(1, 60),
                TIGIT = c(0, 10),
                IgG1 = c(1,10),
                IgG2a = c(1,10),
                IgG2b = c(1,10))



# Build simplex

proteins = rownames(adt_raw)

bg_sim = rep(0, length(proteins))
counts = as.matrix(counts(adt_raw[, colnames(adt_raw) %in% cells0]))

for (i in 1:length(proteins)) {
  
  protein = proteins[i]
  
  value = counts[protein,]
  bg_sim[i] = sum(value[value >= bg_limit[[i]][1] & 
                          value <= bg_limit[[i]][2]])
  
}

names(bg_sim) = proteins
bg_sim = bg_sim/sum(bg_sim)


df = data.frame(mid = mid_bg_profile,
                bg = bg_sim,
                row.names = rownames(adt_raw))


F1E = ggplot(df, aes(mid, bg)) +
  geom_point(shape = 1, size = 2) +
  annotate("text", x=0.35, y=0.01, label=paste0("R = ", round(cor(df$mid, df$bg),3)), size=4) + 
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  geom_text_repel(label = rownames(df)) +
  ggtitle('Cell Background Profile vs. Spongelet Profile') + 
  labs(x = 'Spongelet Profile', y = 'Cell Background Profile') +
  labs(tag = 'E') + 
  coord_cartesian(xlim = c(0,0.4), ylim = c(0, 0.4)) +
  theme_classic() + 
  theme(axis.text = element_text(size = 12),
        plot.tag = element_text(face = "bold",
                                size = 15))

## Save F1E of PBMC5K, MALT10K for FigS8
save(F1E, file = paste0('FigS8_',datadir_name,'.rdata'))




## ====================
## ====== Figure 1
## ====================
l = list()
l[[1]] = F1A
l[[2]] = F1B
l[[3]] = F1C
l[[4]] = F1D
l[[5]] = F1E

design = "
AABBCCC
DDDDEEE
"

F1 = wrap_plots(l) + 
  plot_layout(design = design)

ggsave("Figure1.pdf",
       F1,
       width = 11*1.2,
       height = 8.5*1.2,
       units = "in")


## ===============================
## =========== Supplementary
## ===============================

S2A = DimPlot(adt_seurat, 
              reduction = "adtUMAP",
              label = TRUE,
              repel = TRUE,
              label.size = 3,
              group.by = "cellType") +
  NoLegend() +
  labs(tag = "A") +
  theme(axis.title = element_blank(),
        plot.title = element_blank(),
        plot.tag = element_text())
  


features = rownames(adt_seurat)
S2B = FeaturePlot(adt_seurat, 
                  features = features,
                  combine = FALSE)

for(i in 1:length(S2B)) {
  S2B[[i]] <- S2B[[i]] + NoLegend() + NoAxes() + theme(plot.title = element_text(size=10))
}

S2B = cowplot::plot_grid(plotlist = S2B) +
  labs(tag = "B")



par(mfrow = c(1,1), mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
df =non_zero_drop[cells,]
ggplot(df, aes(rna_size, prot_size, col = adt_seurat$cellType))+
  geom_point() +
  scale_x_continuous(trans = 'log1p', breaks = 10^c(1:5)) +
  scale_y_continuous(trans = 'log1p', breaks = 10^c(1:5)) +
  xlab('RNA Library Size (log)') +
  ylab('ADT Library Size (log)') +
  theme_classic() +
  theme(legend.title=element_blank())



# DE analysis (ADT)
adt.markers <- FindAllMarkers(adt_seurat, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

adt.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
S2C = DoHeatmap(adt_seurat, 
                features = top10$gene,
                group.by = "cellType",
                size = 2.5,
                raster = F) + 
  NoLegend() +
  labs(tag = "C") +
  theme(plot.margin = margin(r = 30, unit = "pt"),
        plot.tag = element_text(face = "bold",
                                size = 15),
        text = element_text(size = 8))



# DE analysis (RNA)
temp = gex_raw[,cells]
rownames(temp) = rowData(temp)$Symbol

rna_seurat = Seurat::CreateSeuratObject(counts(temp))
rna_seurat = Seurat::NormalizeData(rna_seurat)
rna_seurat =  FindVariableFeatures(rna_seurat, selection.method = "vst", nfeatures = 2000)
rna_seurat = ScaleData(rna_seurat)


rna_seurat[["cellType"]] <- adt_seurat[["cellType"]]
Idents(rna_seurat) = Idents(adt_seurat)

rna.markers <- FindAllMarkers(rna_seurat, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
S2D = DoHeatmap(rna_seurat, 
                features = top10$gene,
                group.by = "cellType",
                size = 2.5,
                raster = F) + 
  NoLegend() +
  labs(tag = "D") +
  theme(plot.margin = margin(r = 30, unit = "pt"),
        plot.tag = element_text(face = "bold",
                                size = 15),
        text = element_text(size = 8))


l = list()
l[[1]] = S2A
l[[2]] = S2B
l[[3]] = S2C
l[[4]] = S2D


design = "
ABB
ABB
CCC
DDD
DDD
"

S2 = wrap_plots(l) + 
  plot_layout(design = design)

ggsave("FigS2.pdf",
       S2,
       width = 8.5,
       height = 11,
       units = "in")




## FigS3
library(ggplot2)
library(patchwork)

load('FigS3_data_pbmc5knextgem.rdata')
S3A = F1A + ggtitle('Dataset: PBMC 5K') +
  scale_color_manual(labels = c('Cluster A','Cluster C','Cluster D'),
                     values = cbPalette[c(2,4,5)])

load('FigS3_data_malt10k.rdata')
S3B = F1A + ggtitle('Dataset: MALT 10K') +
  scale_color_manual(labels = c('Cluster A','Cluster B','Cluster C','Cluster D'),
                     values = cbPalette[2:5])

load('FigS3_data_Golomb.rdata')
S3C = F1A + ggtitle('Dataset: Golomb') +
  scale_color_manual(labels = c('Cluster A','Cluster B','Cluster C','Cluster D'),
                     values = cbPalette[2:5])

l = list()
l[[1]] = S3A
l[[2]] = S3B
l[[3]] = S3C



design = "
AA
BB
CC
"

S3 = wrap_plots(l) + 
  plot_layout(design = design)

ggsave("FigS3.pdf",
       S3,
       width = 8.5,
       height = 11,
       units = "in")




## FigS4: Density plot (all ADTs)

max = ifelse(nrow(adt_raw) <= 31, nrow(adt_raw), 31)
proteins = rownames(adt_raw)[1:max]
p = list()

for (i in 1:length(proteins)) {
  
  protein = proteins[i]
  
  df_bg1 = data.frame(Spongelet=assay(adt_raw[, colnames(adt_raw) %in% cells2_all],'counts')[protein,])
  df_bg1.m = reshape2::melt(df_bg1)
  
  df_bg2 = data.frame(Ambient=assay(adt_raw[, colnames(adt_raw) %in% cells3_top],'counts')[protein,])
  df_bg2.m = reshape2::melt(df_bg2)
  
  df_decon <- data.frame(Cells=assay(adt_raw[, colnames(adt_raw) %in% cells0],'counts')[protein,])
  df_decon.m = reshape2::melt(df_decon)
  
  
  p1 = ggplot(df_decon.m, aes(value, fill = variable)) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(trans=scales::pseudo_log_trans(), breaks=c(1,5,10^(1:4))) +
    ggtitle(protein) +
    scale_fill_manual(values=cbPalette[2]) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  
  
  # add background
  ylimit = layer_scales(p1)$y$get_limits()
  ylimit[2] = 1
  
  p1 = p1 + 
    geom_line(data = df_bg1.m, aes(value, color = variable, fill=NA), 
              stat = 'density',
              adjust = 2,
              size = 0.8) + 
    # coord_cartesian(ylim = ylimit) +
    # scale_color_manual(values=cbPalette[c(4,5)]) + 
    theme(legend.title = element_blank(),
          legend.margin=margin(t = -10))
  
  p1 = p1 + 
    geom_line(data = df_bg2.m, aes(value, color = variable, fill=NA), 
              stat = 'density',
              adjust = 2,
              size = 0.8) + 
    coord_cartesian(ylim = ylimit) +
    scale_color_manual(values=cbPalette[c(5,4)]) + 
    theme(legend.title = element_blank(),
          legend.margin=margin(t = -10))
  
  
  p[[i]] = p1
}

S4 = wrap_plots(p, ncol = 4, guides = "collect")+ 
  plot_layout()&theme(legend.position = "bottom", 
                      plot.margin = unit(c(3,3,2,1), "pt"))


save(S4, file = paste0("FigS4_",datadir_name,".rdata"))





# Use saved RData to plot FigS4

datadir_name = 'data_pbmc10khealthdonor'
load(paste0("FigS4_",datadir_name,".rdata"))


datadir_name = 'data_pbmc5knextgem'
load(paste0("FigS4_",datadir_name,".rdata"))


datadir_name = 'data_malt10k'
load(paste0("FigS4_",datadir_name,".rdata"))

datadir_name = 'data_Golomb'
load(paste0("FigS4_",datadir_name,".rdata"))




ggsave(paste0("FigS4_",datadir_name,".pdf"), 
       S4,
       width = 8.5,
       height = 11,
       units = "in")




library(ggplot2)
library(patchwork)

load('FigS8_data_pbmc5knextgem.rdata')
S8A = F1E + ggtitle('Dataset: PBMC 5K') +
  labs(tag = 'A', x = 'Spongelet Profile') +
  annotate("text", x=0.2, y=0.01, label=paste0("R = ", round(cor(F1E$data$mid, F1E$data$bg),3)), size=4) + 
  coord_cartesian(xlim = c(0,0.25), ylim = c(0, 0.25))

load('FigS8_data_malt10k.rdata')
S8B = F1E + ggtitle('Dataset: MALT 10K') +
  labs(tag = 'B', x = 'Spongelet Profile') +
  annotate("text", x=0.2, y=0.01, label=paste0("R = ", round(cor(F1E$data$mid, F1E$data$bg),3)), size=4) + 
  coord_cartesian(xlim = c(0,0.25), ylim = c(0, 0.25))


l = list()
l[[1]] = S8A
l[[2]] = S8B



design = "
A
B
"

S8 = wrap_plots(l) + 
  plot_layout(design = design)

ggsave("FigS8.pdf",
       S8,
       width = 8.5,
       height = 11,
       units = "in")






# Density plot background peak bounds
proteins = rownames(adt_raw)
p = list()

for (i in 1:length(proteins)) {
  protein = proteins[i]
  
  df_decon <- data.frame(value = assay(adt_raw[, colnames(adt_raw) %in% c(cells0)],'counts')[protein,])
  df_decon.m = reshape2::melt(df_decon)
  
  
  # 
  p1 = ggplot(df_decon.m, aes(value)) +
    geom_density(alpha = 0.5, fill = '#F8766D') +
    geom_vline(xintercept = bg_limit[[i]][1]) +
    annotate("text", x = bg_limit[[i]][1], label = as.character(bg_limit[[i]][1]), y=0, vjust = -1) +
    geom_vline(xintercept = bg_limit[[i]][2]) +
    annotate("text", x = bg_limit[[i]][2], label = as.character(bg_limit[[i]][2]), y=0, vjust = -1) +
    scale_x_continuous(trans=scales::pseudo_log_trans(), breaks=c(1,5,10^(1:4))) +
    ggtitle(protein) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title=element_blank())
  
  p[[i]] = p1
  
}

S14 = wrap_plots(p, ncol = 6, guides = "collect")
ggsave(paste0("FigS14_",datadir_name,".pdf"), 
       S14,
       width = 11,
       height = 8.5,
       units = "in")

