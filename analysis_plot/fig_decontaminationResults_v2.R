library(ggplot2)
library(patchwork)
library(scattermore)
library(Seurat)
library(ggridges)
library(celda)
library(tidyverse)
library(cowplot)
library(ggsci)

par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette = rev(cbPalette)
cbPalette = c("#000000", cbPalette)



file_name = 'shrinkage_data_pbmc10khealthdonor'

load(paste0('/rprojectnb2/camplab/home/yin/poisson/model_poisson/',
            file_name,'_Robj.Rdata'))



## ====== Extract optimized param values
val = out@sim$est

r_est = val$r

delta_est = t(val$delta)
delta_mean_est = val$delta_mean

background_est = val$background
background_mean_est = val$background_mean


# Reconstruct decomposed counts

# p
p_est = dat$p

unscaled_rates = matrix(p_est,
                        nrow = nrow(dat$counts),
                        ncol = ncol(dat$counts))

scaling_factor = matrix(dat$OC,
                        nrow = nrow(dat$counts),
                        ncol = ncol(dat$counts),
                        byrow = T)

empty_rate_est = scaling_factor * unscaled_rates * delta_est * (1 - background_est)


# r
unscaled_rates = t(r_est[dat$cell_type,])

scaling_factor = matrix(dat$OC,
                        nrow = nrow(dat$counts),
                        ncol = ncol(dat$counts),
                        byrow = T)

cell_rate_est = scaling_factor * unscaled_rates * (1-delta_est) * (1 - background_est)

# background
background_est = val$background

scaling_factor = matrix(dat$OC,
                        nrow = nrow(dat$counts),
                        ncol = ncol(dat$counts),
                        byrow = T)

background_rate_est = background_est * scaling_factor

rownames(background_est) = rownames(adt_filtered)
rownames(background_rate_est) = rownames(adt_filtered)




# Decontaminated counts
counts = dat$counts

decontaminated_counts = counts * cell_rate_est/(empty_rate_est + cell_rate_est + background_rate_est)
ambient_counts = counts * empty_rate_est/(empty_rate_est + cell_rate_est + background_rate_est)
background_counts = counts * background_rate_est/(empty_rate_est + cell_rate_est + background_rate_est)



## Dot plot of raw counts vs. native rate

prep = function(df, cell_type){
  avg.exp = aggregate(df, by = list(cell_type), FUN = 'mean')
  avg.exp = reshape2::melt(avg.exp, id.vars = c('Group.1'), value.name = 'avg.exp')
  
  ndrop = aggregate(df, by = list(cell_type), FUN = 'length')
  ndrop = reshape2::melt(ndrop,  id.vars = c('Group.1'), value.name = 'ndrop')
  
  ndrop.exp = aggregate(df, by = list(cell_type), FUN = function(x) sum(x > 1))
  ndrop.exp = reshape2::melt(ndrop.exp, id.vars = c('Group.1'),  value.name = 'ndrop.exp')
  
  avg.exp$ndrop = ndrop$ndrop
  avg.exp$ndrop.exp = ndrop.exp$ndrop.exp
  avg.exp$cluster = as.factor(avg.exp$Group.1)
  
  return(avg.exp)
  
}

# DimPlot(adt_seurat, reduction = "adtUMAP", label = TRUE, group.by = "adtClusterID")

cell_type = Idents(adt_seurat)
levels(cell_type) <- list(CD14Mono = "0",
                          CD4NaiveT = "1",
                          CD4MemoryT = "2",
                          NK =  "3", # CD56 CD16
                          CD8MemoryT = "4",
                          B = "5", # CD19 CD20
                          CD8NaiveT = "6",
                          CD16Mono =  "7",
                          Treg =  "8", # CD4 half CD45RA half CD45RO # RNA DE FOXP3 TIGIT: exhausted Treg
                          gdT = "9", # RNA CD3D TRDC TRGC2
                          Doublets =  "10", # Similar to 8
                          Doublets2 =  "11") # B + NK cells


raw = data.frame(t(counts))
raw_df = prep(raw, cell_type)


decon = data.frame(t(decontaminated_counts))
decon_df = prep(decon, cell_type)


F3A = raw_df %>% 
  mutate(`% Expressing` = (ndrop.exp/ndrop) * 100) %>% 
  filter(avg.exp > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = variable, color = avg.exp, size = `% Expressing`)) + 
  geom_point() +
  scale_color_viridis_c(name = 'Average\ncounts',
                        trans = 'log1p',
                        breaks = 10^c(1:6)) + 
  scale_size(range = c(1,5)) +
  cowplot::theme_cowplot() + 
  ggtitle('Marker Expression in Original Counts') + 
  labs(x = "", y = "", size = "%\nExpressing", tag = "A") +
  theme(axis.line  = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        plot.tag = element_text(face = "bold",
                                size = 15)) +
  guides(size = guide_legend(order = 1))


F3B = decon_df %>% 
  mutate(`% Expressing` = (ndrop.exp/ndrop) * 100) %>% 
  filter(avg.exp > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = variable, color = avg.exp, size = `% Expressing`)) + 
  geom_point() +
  scale_color_viridis_c(name = 'Average\ncounts',
                        trans = 'log1p',
                        breaks = 10^c(1:6)) + 
  scale_size(range = c(1,5)) +
  cowplot::theme_cowplot() + 
  ggtitle('Marker Expression in Decontaminated Counts') + 
  labs(x = "", y = "", size = "%\nExpressing", tag = "B") +
  theme(axis.line  = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        plot.tag = element_text(face = "bold",
                                size = 15)) +
  guides(size = guide_legend(order = 1))



## Background density (Zoom in)

proteins = c("CD3", "CD16")
clusters = c("CD8NaiveT", "CD16Mono", "Treg", "NK")

counts = counts(adt_filtered)
counts_norm = sweep(counts, 2, colSums2(counts), '/')


p = list()

for(i in 1:length(proteins)){
  
  protein = proteins[i]
  
  
  df_ridge <- data.frame(raw_counts_scaled = counts_norm[protein,],
                         background_estimated = background_est[protein,],
                         cell_type = cell_type)
  df_ridge = df_ridge[df_ridge$cell_type %in% clusters, ]
  df_ridge.m = reshape2::melt(df_ridge)
  
  
  
  p1 = ggplot(df_ridge.m, aes(x = value, y = cell_type, fill = variable)) +
    geom_density_ridges(scale = 1, alpha = 0.5) +
    scale_fill_manual(values = cbPalette[4:5], labels = c('Normalized Counts', 'Est. Background')) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic() +
    ggtitle(paste0(protein)) + 
    labs(x = "", y = "") +
    coord_cartesian(xlim = c(0,0.4)) +
    theme(legend.position="bottom",
          legend.title = element_blank(),
          axis.text.y = element_text(angle = 90,hjust = 0.1,
                                     size = 10),
          legend.margin=margin(t = -10))
  
  
  
  
  p[[i]] = p1
  
}

p[[1]] = p[[1]] + labs(tag = "C") 

F3C = wrap_plots(p, nrow = 1, guides = "collect") + 
  plot_layout()&theme(legend.position = "bottom", 
                      legend.text = element_text(size = 10),
                      legend.key.size = unit(0.5, 'cm'),
                      plot.margin = unit(c(3,3,2,1), "pt"),
                      axis.text.x = element_text(angle = 90),
                      axis.title.y = element_text(margin = margin(r = -10)),
                      plot.tag = element_text(face = "bold",
                                              size = 15))




## Decontamination barplot
cellTypeMappings <- list(MonoCD14 = 0,
                         CD4T = c(1,2),
                         NK = 3,
                         CD8T = c(4,6),
                         B = 5,
                         MonoCD16 = 7)

sce = SingleCellExperiment(list(counts = dat$counts,
                                decon = decontaminated_counts))
colData(sce)$decontX_clusters = adt_seurat$adtClusterID


markers <- list(CD14 = c("CD14"),
                CD4 = c("CD4"),
                CD56 = c("CD56"),
                CD8 = c("CD8a"),
                CD19 = c("CD19"),
                CD16 = c("CD16")

)


F3D_temp = plotDecontXMarkerPercentage(sce,
                            markers = markers,
                            groupClusters = cellTypeMappings,
                            assayName = c("counts", "decon"),
                            ncol = 2,
                            labelSize = 2.5) 

data_temp = F3D_temp$data

F3D = ggplot(data_temp, aes(markerLabels, percent, fill = assay)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(x = markerLabels, y = percent + 2.5, label = percent),
            position = ggplot2::position_dodge2(width = 0.9, preserve = "single"),
            size = 2.5) +
  facet_wrap(vars(cellTypeLabels), ncol = 2) +
  labs(y = 'Percentage of cells expressing markers') +
  scale_fill_discrete(labels = c('Original', 'Decontaminated')) +
  labs(tag = "D") +
  theme(legend.position = "bottom",
        legend.margin=margin(t = 0),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1,
                                   size = 10),
        plot.tag = element_text(face = "bold",
                                size = 15),
        panel.background = element_rect(
          fill = "white",
          color = "grey"
        ),
        panel.grid = element_line("grey"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 10))
  


## ========= Percentage of native, ambient, background

df = data.frame(Native = colSums(decontaminated_counts),
                Ambient = colSums(ambient_counts),
                Background = colSums(background_counts))

df = sweep(df, 1, rowSums(df), '/')

df.m = reshape2::melt(df)

F3E = ggplot(df.m, aes(x = variable, y = value, fill = variable)) +
  geom_jitter(size = 0.1, color = 'gray') +
  geom_violin(alpha = 0.8) +
  # scale_y_continuous(trans = "log1p", breaks = 10^c(1:6))+ 
  scale_fill_brewer(palette="RdBu") +
  labs(y = 'Proportion', tag = "E") +
  ggtitle('Proportion of Counts from Different Sources') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        plot.tag = element_text(face = "bold",
                                size = 15))





## ===== UMAP before and after decontamination

# old umap
adt_seurat[["adtClusterID_name"]] <- cell_type

F3F1 = DimPlot(adt_seurat, 
             reduction = "adtUMAP", 
             label = TRUE,
             label.size = 2.5,
             repel = TRUE,
             group.by = "adtClusterID_name") +
  labs(tag = "F") +
  ggtitle("Before Decontamination") +
  theme(plot.title = element_text(face = "plain",
                                  size = 12),
        plot.tag = element_text(face = "bold",
                                size = 15),
        axis.title = element_text(size = 10))


adt_seurat_decon = Seurat::CreateSeuratObject(decontaminated_counts,
                                              assay = 'ADT')
adt_seurat_decon <- Seurat::NormalizeData(adt_seurat_decon, 
                                    normalization.method = "CLR",
                                    margin = 2)
adt_seurat_decon <- ScaleData(adt_seurat_decon, assay = "ADT", verbose = FALSE)
adt_seurat_decon <- RunPCA(adt_seurat_decon,
                           assay = "ADT",
                           features = rownames(adt_seurat_decon),
                           reduction.name = "pca_adt",
                           reduction.key = "pcaadt_", 
                           verbose = FALSE,
                           npcs = 10)

adt_seurat_decon <- RunUMAP(adt_seurat_decon,
                            dims = 1:10,
                            assay = "ADT",
                            reduction = "pca_adt",
                            reduction.name = "adtUMAP",
                            verbose = FALSE)

adt_seurat_decon[["adtClusterID"]] <- cell_type # adt_seurat$adtClusterID

F3F2 = DimPlot(adt_seurat_decon,
             reduction = "adtUMAP",
             label = TRUE,
             label.size = 2.5, 
             repel = TRUE,
             group.by = "adtClusterID") +
  ggtitle('After Decontamination') +
  theme(plot.title = element_text(face = "plain",
                                  size = 12),
        axis.title = element_text(size = 10))

F3F = F3F1 + F3F2+
  plot_layout(nrow=1, widths=c(1,1), guides = 'collect')&theme(legend.position = "none")




## ====================
## ====== Figure 3
## ====================
l = list()
l[[1]] = F3A
l[[2]] = F3B
l[[3]] = F3C
l[[4]] = F3D
l[[5]] = F3E
l[[6]] = F3F

design = "
AAAAAABBBBBB
AAAAAABBBBBB
CCCCDDDEEEEE
CCCCDDDFFFFF
CCCCDDDFFFFF
"

F3 = wrap_plots(l) + 
  plot_layout(design = design)&theme(legend.margin=margin(t = -30, b = 20))

ggsave("Figure3.pdf",
       F3,
       width = 11*1.3,
       height = 8.5*1.3,
       units = "in")








## ===============================
## =========== Supplementary
## ===============================


## Density plot
proteins = rownames(adt_filtered)

p = list()

for (i in 1:length(proteins)) {
  
  protein = proteins[i]
  
  
  df_poi <- data.frame(counts=assay(adt_filtered,'counts')[protein,],
                       decon=decontaminated_counts[protein,])
  df_poi.m = reshape2::melt(df_poi)
  
  
  
  # Poisson plot
  p1 = ggplot(df_poi.m, aes(value, fill = variable)) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(trans=scales::pseudo_log_trans(), breaks=c(1,5,10^(1:4))) +
    # scale_fill_discrete() +
    scale_fill_npg(alpha = 0.7, labels = c('Original', 'Decontaminated')) + 
    ggtitle(protein) +
    labs(x = "", fill ="PBMC10k") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9),
          # legend.title = element_blank(),
          legend.margin=margin(t = -10))
  
  ylimit = layer_scales(p1)$y$get_limits()
  ylimit[2] = min(ylimit[2], 5)
  
  p[[i]] = p1 + coord_cartesian(ylim = ylimit)
  
  
}

S13A = wrap_plots(p, ncol = 6, guides = "collect") +
  plot_layout()&theme(legend.position = "bottom", 
                      plot.margin = unit(c(3,3,2,1), "pt"))



# Violin plot of B cell after decomtam
df = data.frame(ori = counts(adt_filtered)['CD14',],
                decom = decontaminated_counts['CD14',],
                cell_type = cell_type)

df.m = reshape2::melt(df)

S13B = ggplot(df.m, aes(x = cell_type, y = value, fill = variable)) + 
  geom_violin(scale = 'width') + 
  scale_fill_discrete(labels = c('Original', 'Decontaminated')) + 
  ylim(0,250) +
  cowplot::theme_cowplot() + 
  ggtitle('CD14 Counts by Cell Types') + 
  labs(x = "", y = "Counts", tag = "B") +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9),
        legend.title = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 11),
        legend.margin=margin(t = -30, l=150), 
        plot.tag = element_text(face = "bold",
                                size = 15))





l = list()
l[[1]] = S13A
l[[2]] = S13B



design = "
A
A
B
"

S13 = wrap_plots(l) + plot_layout(design = design)

ggsave("FigS13.pdf",
       S13,
       width = 8.5,
       height = 11,
       units = "in")
