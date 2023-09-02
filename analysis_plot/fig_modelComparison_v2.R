library(Seurat)
library(ggplot2)
library(patchwork)
library(scater)
library(cluster)
library(dplyr)
library(ggsci)
library(stringr)
library(ggridges)


cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette = rev(cbPalette)
cbPalette = c("#000000", cbPalette)



## =========== Import all datasets

all_datasets = list()


## 10x datasets
datasets = c('data_pbmc5knextgem', 
             'data_pbmc10khealthdonor',
             'data_malt10k',
             'young_aged')

for (i in 1:length(datasets)) {
  
  dataname = datasets[i]
  
  file_poi = paste0('shrinkage_',dataname,'_poisson.Rdata')
  file_nb = paste0('shrinkage_',dataname,'_nb.Rdata')
  file_dsb = paste0(dataname,'_dsb.Rdata')
  file_totalvi = paste0('totalvi_denoised_',dataname,'.csv')
  file_scar = paste0('scAR_denoised_',dataname,'.csv')
  
  temp = list()
  
  # Poisson
  
  # 10x dir
  poi_dir = paste0('/rprojectnb2/camplab/home/yin/poisson/model_poisson/', file_poi)
  if (file.exists(poi_dir)) {
    load(poi_dir)
  } else {
    # Zhang dir
    poi_dir = paste0('/restricted/projectnb/camplab/home/yin/poisson/analysis_zhang/', file_poi)
    if (file.exists(poi_dir)) {
      load(poi_dir)
    }
    
  }
  
  
  temp[['gex_filtered']] = gex_filtered
  temp[['adt_filtered']] = adt_filtered
  temp[['adt_seurat']] = adt_seurat
  temp[['poi']] = decontaminated_counts
  temp[['decontX']] = assay(adt_filtered, 'decontXcounts')
  
  
  ## nb
  load(paste0('/restricted/projectnb/camplab/home/yin/poisson/submission_nar_NB_vs_poi/', file_nb))
  temp[['gex_filtered_nb']] = gex_filtered
  temp[['adt_filtered_nb']] = adt_filtered
  temp[['adt_seurat_nb']] = adt_seurat
  temp[['nb']] = decontaminated_counts
  
  
  ## dsb
  load(paste0('/rprojectnb2/camplab/home/yin/poisson/method_dsb/', file_dsb))
  cell.adt.mtx = cell.adt.mtx[!grepl('IgG', rownames(cell.adt.mtx)),]
  cell.adt.mtx = cell.adt.mtx[!grepl('HTO', rownames(cell.adt.mtx)),]
  # rownames(cell.adt.mtx) = str_remove(rownames(cell.adt.mtx), "_TotalSeqB")
  # rownames(cell.adt.mtx) = str_remove(rownames(cell.adt.mtx), "CITE_")
  cell.dsb.norm = cell.dsb.norm[!grepl('IgG', rownames(cell.dsb.norm)),]
  cell.dsb.norm = cell.dsb.norm[!grepl('HTO', rownames(cell.dsb.norm)),]
  # rownames(cell.dsb.norm) = str_remove(rownames(cell.dsb.norm), "_TotalSeqB")
  # rownames(cell.dsb.norm) = str_remove(rownames(cell.dsb.norm), "CITE_")
  
  rownames(cell.adt.mtx) = rownames(adt_filtered)
  rownames(cell.dsb.norm) = rownames(adt_filtered)
  
  temp[['raw']] = cell.adt.mtx
  temp[['dsb']] = cell.dsb.norm
  
  
  
  
  
  ## totalVI
  totalvi_denoised = read.csv(paste0('/rprojectnb2/camplab/home/yin/poisson/method_totalVI/',file_totalvi),
                              header = T,
                              row.names = 1)
  totalvi_denoised = totalvi_denoised[,!grepl('IgG', colnames(totalvi_denoised))]
  totalvi_denoised = totalvi_denoised[,!grepl('HTO', colnames(totalvi_denoised))]
  # colnames(totalvi_denoised) = str_remove(colnames(totalvi_denoised), "_TotalSeqB")
  # colnames(totalvi_denoised) = str_remove(colnames(totalvi_denoised), "CITE_")
  # totalvi_denoised = totalvi_denoised[,!(colnames(totalvi_denoised) %in% c("CD335", "CXCR4", "KLRG1", "CITE.CD49d"))]
  totalvi_denoised =t(as.matrix(totalvi_denoised))
  rownames(totalvi_denoised) = rownames(adt_filtered)
  
  
  temp[['totalvi']] = totalvi_denoised
  
  
  
  ## scAR
  scar_denoised = read.csv(paste0('/rprojectnb2/camplab/home/yin/poisson/method_scAR/',file_scar),
                           header = T,
                           row.names = 1)
  
  scar_denoised = scar_denoised[,!grepl('IgG', colnames(scar_denoised))]
  scar_denoised = scar_denoised[,!grepl('HTO', colnames(scar_denoised))]
  # colnames(scar_denoised) = str_remove(colnames(scar_denoised), "_TotalSeqB")
  # colnames(scar_denoised) = str_remove(colnames(scar_denoised), "CITE_")
  scar_denoised =t(as.matrix(scar_denoised))
  rownames(scar_denoised) = rownames(adt_filtered)
  
  
  temp[['scar']] = scar_denoised
  
  
  all_datasets[[dataname]] = temp
  
  
}
                    








## ======= Annotate cluster
# PBMC10k
# features = c("CD3", "CD4", "CD8a","CD19","CD14", "CD16", "CD56","CD127", "CD45RA", "CD45RO", "PD-1", "TIGIT")
# PBMC5k
# features = c("CD3", "CD4", "CD8a","CD19","CD20","CD14", "CD16", "CD45RA", "CD45RO", "CD56", "TIGIT")
# MALT10k
# features = rownames(adt_seurat)
# Zhang
features = c("Tmem119", "Mrc1", # microglia: Tmem119+,Mrc1–
  "Cd79a", "Ly6d", "Iglc2", # Bcell
  "Thy1", "Itga2", "Klrb1" # T: Thy1+/Itga2+/Klrb1-
  # "Plac8", "Fn1", "S100a4", # Ly6cHigh
  # "H2-Aa", "H2-Eb1", "Cd74" # Ly6CLow
  # " S100a9", "S100a8", "Retnlg", "Lcn2" # Ly6cHigh Ly6CLow
) 


features = c('CD3',
             'CD8a',
             'CD11b', # microglia: CD45Low,CD11b+,CD38Low,MHC-IILow
             'CD38', # BAM
             'CD45', # peripherally derived immune cells
             'Ly6C',
             'Ly6G',
             'CD45R-B220', # B cell
             'NK1.1' # NK cell
)

# Feature plot
# umap = adt_seurat@reductions$adtUMAP@cell.embeddings

p1 = DimPlot(adt_seurat, reduction = "adtUMAP", label = TRUE, group.by = "adtClusterID")
# p2 = FeaturePlot(adt_seurat, features = features, combine = F)
# 
# for(i in 1:length(p2)) {
#   p2[[i]] <- p2[[i]] + NoLegend() + NoAxes()
# }
# 
# p2 = cowplot::plot_grid(plotlist = p2)
# 
# p1 + p2 + plot_layout(nrow=1,widths=c(1,2))

p2 = FeaturePlot(adt_seurat, features = features)

p1 + p2 + plot_layout(nrow=1,widths=c(1,3))





## DE analysis on ADT space

adt.markers <- FindAllMarkers(adt_seurat, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)


adt.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(adt_seurat, features = top10$gene) + NoLegend()




## DE on gex_filtered to identify clusters

rna_seurat = Seurat::CreateSeuratObject(counts(gex_filtered))
rna_seurat <- Seurat::NormalizeData(rna_seurat)
rna_seurat =  FindVariableFeatures(rna_seurat, selection.method = "vst", nfeatures = 2000)
rna_seurat = ScaleData(rna_seurat)
# rna_seurat <- RunPCA(rna_seurat, features = VariableFeatures(rna_seurat))
# rna_seurat <- FindNeighbors(rna_seurat, dims = 1:50)
# rna_seurat <- FindClusters(rna_seurat, resolution = 0.8)
# rna_seurat <- RunUMAP(rna_seurat, dims = 1:50)

# Add ADT UMAP
rna_seurat@reductions$adtUMAP = adt_seurat@reductions$adtUMAP
# Add ADT cluster labels
rna_seurat[["adtClusterID"]] <- Idents(adt_seurat)
Idents(rna_seurat) = Idents(adt_seurat)

rna.markers <- FindAllMarkers(rna_seurat, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)


rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(rna_seurat, features = top10$gene) + NoLegend()



# For PBMC10k
# gamma delta T cell: CD3D TRDC TRGC2
VlnPlot(rna_seurat, features = c('FCER1A'))

cluster.markers <- FindMarkers(rna_seurat,
                               ident.1 = 9,
                               ident.2 = c(2, 4),
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

cluster.markers <- FindMarkers(rna_seurat,
                               ident.1 = 10,
                               ident.2 = 8,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)



cluster.markers %>%
  top_n(n = 10, wt = avg_log2FC) -> cluster.top10

# Cluster 10
VlnPlot(rna_seurat, features = c('CD3D'))
VlnPlot(rna_seurat, features = c('CD4'))
VlnPlot(rna_seurat, features = c('NKG7'))
VlnPlot(rna_seurat, features = c('CD14'))
VlnPlot(rna_seurat, features = c('LYZ'))
VlnPlot(rna_seurat, features = c('CD163'))
VlnPlot(rna_seurat, features = c('MS4A1'))




# For MALT10k
VlnPlot(rna_seurat, features = c('FCER1A'))
VlnPlot(rna_seurat, features = c('IRF7'))


cluster.markers <- FindMarkers(rna_seurat,
                                ident.1 = 7,
                                ident.2 = c(2),
                                min.pct = 0.25,
                                logfc.threshold = 0.25)

cluster.markers %>%
  top_n(n = 10, wt = avg_log2FC) -> cluster.top10






p1 = DimPlot(rna_seurat, reduction = 'adtUMAP', label = TRUE, group.by = "adtClusterID")
# malt 10k
features = c('FOXP3', 'RTKN2', 'TIGIT', 'CXCL13', 'PDCD1')
p2 = FeaturePlot(rna_seurat, reduction = 'adtUMAP', features = features)

p1 + p2 + plot_layout(nrow=1,widths=c(1,2))




## ======== Annotate cluster in Seurat obj

# Label for pbmc5knextgem
cellTypeMappings <- list(MonoCD14 = "0",
                         CD4TNaive = "1",
                         CD4TMemory = "2",
                         CD8TMemory =  "3",
                         NK = "4",
                         B = "5")

# cell_type mapping
cell_type = all_datasets$data_pbmc5knextgem$adt_seurat$adtClusterID
levels(cell_type) = cellTypeMappings
all_datasets$data_pbmc5knextgem$adt_seurat[["cellType"]] = cell_type




# Label for pbmc10khealthdonor
cellTypeMappings <- list(CD14Mono = "0",
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

# cell_type mapping
cell_type = all_datasets$data_pbmc10khealthdonor$adt_seurat$adtClusterID
levels(cell_type) = cellTypeMappings
all_datasets$data_pbmc10khealthdonor$adt_seurat[["cellType"]] = cell_type




# Label for malt10k
cellTypeMappings <- list(Bcell_IGHD = "0", # IGHD, FCER2 (violin plot)
                         Bcell_IGHM = "1", # IGHM, IGHA1
                         CD4TMemory = "2",
                         CD8TMemory = "3",
                         Treg = "4", # 'FOXP3', 'RTKN2', 'TIGIT'
                         Tfh = "5", # 'TIGIT', 'CXCL13', 'PDCD1'
                         doublets = "6", #?
                         uk2 = "7",
                         DC = "8",
                         uk3 = "9")

# cell_type mapping
cell_type = all_datasets$data_malt10k$adt_seurat$adtClusterID
levels(cell_type) = cellTypeMappings
all_datasets$data_malt10k$adt_seurat[["cellType"]] = cell_type




# Label for Zhang
cellTypeMappings <- list(Microglia = "0",
                         CD8T = "2",
                         B = "5",
                         CD4T = "8")


lognom = function(mat){
  mat_norm = sweep(mat, 2, colSums(mat), '/')
  mat_norm = log1p(10000*mat_norm)
}


clrnorm = function(mat){
  temp = Seurat::CreateSeuratObject(mat,
                                    assay = 'ADT')
  temp <- Seurat::NormalizeData(temp, 
                                normalization.method = "CLR",
                                margin = 2)
  temp[['ADT']]@data
}


for (i in 1:length(all_datasets)) {
  all_datasets[[i]]$raw_norm = clrnorm(counts(all_datasets[[i]]$adt_filtered))
  all_datasets[[i]]$scar_denoised_norm = clrnorm(all_datasets[[i]]$scar)
  all_datasets[[i]]$totalvi_denoised_norm = clrnorm(all_datasets[[i]]$totalvi)
  all_datasets[[i]]$poi_norm = clrnorm(all_datasets[[i]]$poi)
}




## ==================== Silhouette width

# Mean Silhouette Width for each clusterID
sil = function(id, matrix){
  sil_score = silhouette(id, dist(t(matrix)))
  
  mean_sil_wid = data.frame(sil_score) %>%
    group_by(cluster) %>%
    summarise(mean = mean(sil_width))
}


## mean silhouette width difference

mean_sil_wid_diff = data.frame()

for (i in 1:length(datasets)){
  
  dataname = datasets[i]
  dat = all_datasets[[dataname]]
  
  clusterIDs = as.integer(dat$adt_seurat$adtClusterID)
  # clusterIDs = as.integer(dat$rna_seurat$seurat_clusters)
  
  sil_old = sil(clusterIDs, dat$raw_norm)
  sil_poi = sil(clusterIDs, dat$poi_norm)
  sil_dsb = sil(clusterIDs, dat$dsb)
  sil_totalvi = sil(clusterIDs, dat$totalvi_denoised_norm)
  sil_scar = sil(clusterIDs, dat$scar_denoised_norm)
  
  
  temp = data.frame(cluster = factor(sil_old$cluster),
                    dataset = dataname, 
                    Poisson = sil_poi$mean - sil_old$mean,
                    dsb = sil_dsb$mean - sil_old$mean,
                    totalVI = sil_totalvi$mean - sil_old$mean,
                    scAR = sil_scar$mean - sil_old$mean)
  
  
  mean_sil_wid_diff = rbind(mean_sil_wid_diff, temp)
  print(paste0(dataname, " done!"))
  
  
}


mean_sil_wid_diff.m = reshape2::melt(mean_sil_wid_diff)
mean_sil_wid_diff.m$dataset = factor(mean_sil_wid_diff.m$dataset,
                                     levels = c("data_pbmc10khealthdonor",
                                                "data_pbmc5knextgem",
                                                "data_malt10k",
                                                "young_aged"))
levels(mean_sil_wid_diff.m$dataset) = c('PBMC10k', 'PBMC5k', 'MALT10k', 'Golomb')


F4A = ggplot(mean_sil_wid_diff.m, aes(x=dataset, y = value, fill = variable)) +
  geom_boxplot() +
  # geom_jitter(aes(color = cluster, fill = variable)) +
  geom_point(position = position_jitterdodge(0.2)) + 
  scale_fill_discrete(labels = c('DecontPro', 'dsb', 'totalVI', 'scAR')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = 20),
        plot.margin = margin(0,2,0,0),
        legend.margin = margin(-40,0,0,0),
        legend.position = "bottom",
        plot.tag = element_text(face = "bold",
                                size = 15)) + #angle = 45, vjust = 1, hjust = 0.9)))
  labs(y = 'Difference in Mean Silhouette Width within Clusters',
       fill = "",
       tag = "A") +
  ggtitle('Silhouette Width')







## ========== Scoring decontamination results
## Percentage of droplets in defined clusters expressing markers

# pbmc10khealthdonor
# 
cluster_markers_pbmc10k = list(CD14_Monocytes = list(0, c('CD14')),
                               CD4_naive_Tcell = list(1, c('CD4', 'CD45RA')),
                               CD4_memory_Tcell = list(2, c('CD4', 'CD45RO')),
                               Nk =  list(3, c('CD16', 'CD56')),
                               CD8_activated_Tcell = list(4, c('CD8a', 'CD45RO')),
                               Bcell = list(5, c('CD19')),
                               CD8_naive_Tcell = list(6, c('CD8a', 'CD45RA')),
                               CD16_Monocytes = list(7, c('CD16')))

# markers that should NOT be expressed
cluster_non_markers_pbmc10k = list(CD4_naive_Tcell = list(1, c('CD8a', 'CD45RO')),
                                   CD4_memory_Tcell = list(2, c('CD8a', 'CD45RA')),
                                   CD8_activated_Tcell = list(4, c('CD4', 'CD45RA')),
                                   Bcell = list(5, c('CD3')),
                                   CD8_naive_Tcell = list(6, c('CD4', 'CD45RO')))






# pbmc5k
# 
cluster_markers_pbmc5k = list(CD14_Monocytes = list(0, c('CD14')),
                              CD4_naive_Tcell = list(1, c('CD4', 'CD45RA')),
                              CD4_memory_Tcell = list(2, c('CD4', 'CD45RO')),
                              CD8_activated_Tcell = list(3, c('CD8a', 'CD45RO')),
                              Nk =  list(4, c('CD16', 'CD56')), 
                              Bcell = list(5, c('CD19')))

# markers that should NOT be expressed
cluster_non_markers_pbmc5k = list(CD4_naive_Tcell = list(1, c('CD8a', 'CD45RO')),
                                  CD4_memory_Tcell = list(2, c('CD8a', 'CD45RA')),
                                  CD8_activated_Tcell = list(3, c('CD4', 'CD45RA')),
                                  Bcell = list(5, c('CD3')))





# malt10k
# 
cluster_markers_malt10k = list(Bcell1 = list(0, c('CD19')),
                               Bcell2 = list(1, c('CD19')),
                               CD4_memory_Tcell = list(2, c('CD4', 'CD45RO')),
                               CD8_activated_Tcell = list(3, c('CD8a', 'CD45RO')))

# markers that should NOT be expressed
cluster_non_markers_malt10k = list(Bcell1 = list(0, c('CD3')),
                                   Bcell2 = list(1, c('CD3')),
                                   CD4_memory_Tcell = list(2, c('CD8a', 'CD45RA')),
                                   CD8_activated_Tcell = list(3, c('CD4', 'CD45RA')))



# zhang
# 
cluster_markers_zhang = list(Microglia = list(0, c('CD11b')),
                             CD8_Tcell = list(2, c('CD3', 'CD8a')),
                             Bcell = list(5, c('CD45R-B220')))

# markers that should NOT be expressed
cluster_non_markers_zhang = list(Microglia = list(0, c('CD3', 'CD8a')),
                             CD8_Tcell = list(2, c('CD45R-B220')),
                             Bcell = list(5, c('CD3', 'CD8a')))




# Helper function
# Percentage of drops expressing markers
percentageExpressingMarker = function(mat,
                                      cell_type,
                                      clusterID,
                                      markers,
                                      threshold){
  # cell_type: vector w/ length = ncol of mat 
  # clusterID: number
  # markers: vector of marker names
  # threshold: above which to be deemed being expressed
  
  
  temp = mat[markers, cell_type == clusterID]
  temp = as.matrix(temp) > threshold
  
  # Drops that express all markers
  if (dim(temp)[2] == 1) {
    expressingAll = temp
  } else {
    expressingAll = apply(temp, 2, all)
  }
  
  out = c(sum(expressingAll), length(expressingAll))
  return(out)
}






# Calc score for each datasets, iterating algorithms
score_datasets = function(dat_list, 
                          nPosDrops,
                          nDrops,
                          cluster_markers){
  
  
  for(i in 1:length(cluster_markers)){
    
    clusterID = as.character(cluster_markers[[i]][[1]])
    markers = cluster_markers[[i]][[2]]
    
    out = percentageExpressingMarker(dat_list$raw,
                                     dat_list$adt_seurat$adtClusterID,
                                     clusterID,
                                     markers,
                                     1)
    
    nPosDrops["raw"] = nPosDrops["raw"] + out[1]
    nDrops["raw"] = nDrops["raw"] + out[2]
    
    
    out = percentageExpressingMarker(dat_list$decontX,
                                     dat_list$adt_seurat$adtClusterID,
                                     clusterID,
                                     markers,
                                     1)
    
    nPosDrops["decontX"] = nPosDrops["decontX"] + out[1]
    nDrops["decontX"] = nDrops["decontX"] + out[2]
    
    
    out = percentageExpressingMarker(dat_list$nb,
                                     dat_list$adt_seurat$adtClusterID,
                                     clusterID,
                                     markers,
                                     1)
    
    nPosDrops["nb"] = nPosDrops["nb"] + out[1]
    nDrops["nb"] = nDrops["nb"] + out[2]
    
    
    
    out = percentageExpressingMarker(dat_list$poi,
                                     dat_list$adt_seurat$adtClusterID,
                                     clusterID,
                                     markers,
                                     1)
    
    nPosDrops["poi"] = nPosDrops["poi"] + out[1]
    nDrops["poi"] = nDrops["poi"] + out[2]
    
    
    
    out = percentageExpressingMarker(dat_list$dsb,
                                     dat_list$adt_seurat$adtClusterID,
                                     clusterID,
                                     markers,
                                     1)
    
    nPosDrops["dsb"] = nPosDrops["dsb"] + out[1]
    nDrops["dsb"] = nDrops["dsb"] + out[2]
    
    
    
    out = percentageExpressingMarker(dat_list$dsb,
                                     dat_list$adt_seurat$adtClusterID,
                                     clusterID,
                                     markers,
                                     5)
    
    nPosDrops["dsb_low"] = nPosDrops["dsb_low"] + out[1]
    nDrops["dsb_low"] = nDrops["dsb_low"] + out[2]
    
    
    
    out =percentageExpressingMarker(dat_list$totalvi,
                                    dat_list$adt_seurat$adtClusterID,
                                    clusterID,
                                    markers,
                                    1)
    
    nPosDrops["totalvi"] = nPosDrops["totalvi"] + out[1]
    nDrops["totalvi"] = nDrops["totalvi"] + out[2]
    
    
    
    out =percentageExpressingMarker(dat_list$scar,
                                    dat_list$adt_seurat$adtClusterID,
                                    clusterID,
                                    markers,
                                    1)
    
    nPosDrops["scar"] = nPosDrops["scar"] + out[1]
    nDrops["scar"] = nDrops["scar"] + out[2]
    
  }
  
  
  return(list(nPosDrops, nDrops))
  
}





# Scores for PBMC10K
nPosDrops_pbmc10k = c(raw=0, decontX=0, nb=0, poi=0, dsb = 0, dsb_low = 0, totalvi = 0, scar = 0)
nDrops_pbmc10k = c(raw = 0, decontX=0, nb=0, poi=0, dsb = 0, dsb_low = 0, totalvi = 0, scar = 0)

out = score_datasets(all_datasets$data_pbmc10khealthdonor,
                     nPosDrops_pbmc10k,
                     nDrops_pbmc10k,
                     cluster_markers_pbmc10k)

score_pbmc10k = out[[1]]/out[[2]]

out = score_datasets(all_datasets$data_pbmc10khealthdonor,
                     nPosDrops_pbmc10k,
                     nDrops_pbmc10k,
                     cluster_non_markers_pbmc10k)

non_score_pbmc10k = out[[1]]/out[[2]]

# Overall score as positive score minus negative score
overall_score_pbmc10k = score_pbmc10k - non_score_pbmc10k




# Scores for PBMC5K
nPosDrops_pbmc5k = c(raw = 0, decontX=0, nb=0, poi=0, dsb = 0, dsb_low = 0, totalvi = 0, scar = 0)
nDrops_pbmc5k = c(raw = 0, decontX=0, nb=0, poi=0, dsb = 0, dsb_low = 0, totalvi = 0, scar = 0)

out = score_datasets(all_datasets$data_pbmc5knextgem,
                     nPosDrops_pbmc5k,
                     nDrops_pbmc5k,
                     cluster_markers_pbmc5k)

score_pbmc5k = out[[1]]/out[[2]]

out = score_datasets(all_datasets$data_pbmc5knextgem,
                     nPosDrops_pbmc5k,
                     nDrops_pbmc5k,
                     cluster_non_markers_pbmc5k)

non_score_pbmc5k = out[[1]]/out[[2]]

# Overall score as positive score minus negative score
overall_score_pbmc5k = score_pbmc5k - non_score_pbmc5k





# Scores for MALT10K
nPosDrops_malt10k = c(raw = 0, decontX=0, nb=0, poi=0, dsb = 0, dsb_low = 0, totalvi = 0, scar = 0)
nDrops_malt10k = c(raw = 0, decontX=0, nb=0, poi=0, dsb = 0, dsb_low = 0, totalvi = 0, scar = 0)

out = score_datasets(all_datasets$data_malt10k,
                     nPosDrops_malt10k,
                     nDrops_malt10k,
                     cluster_markers_malt10k)

score_malt10k = out[[1]]/out[[2]]

out = score_datasets(all_datasets$data_malt10k,
                     nPosDrops_malt10k,
                     nDrops_malt10k,
                     cluster_non_markers_malt10k)

non_score_malt10k = out[[1]]/out[[2]]

# Overall score as positive score minus negative score
overall_score_malt10k = score_malt10k - non_score_malt10k




# Scores for Zhang
nPosDrops_zhang = c(raw = 0, decontX=0, nb=0, poi=0, dsb = 0, dsb_low = 0, totalvi = 0, scar = 0)
nDrops_zhang = c(raw = 0, decontX=0, nb=0, poi=0, dsb = 0, dsb_low = 0, totalvi = 0, scar = 0)

out = score_datasets(all_datasets$young_aged,
                     nPosDrops_zhang,
                     nDrops_zhang,
                     cluster_markers_zhang)

score_zhang = out[[1]]/out[[2]]

out = score_datasets(all_datasets$young_aged,
                     nPosDrops_zhang,
                     nDrops_zhang,
                     cluster_non_markers_zhang)

non_score_zhang = out[[1]]/out[[2]]


# Overall score as positive score minus negative score
overall_score_zhang = score_zhang - non_score_zhang






## Plot scores
df = data.frame(PBMC10k = c(score_pbmc10k, non_score_pbmc10k, overall_score_pbmc10k),
                PBMC5k = c(score_pbmc5k, non_score_pbmc5k, overall_score_pbmc5k),
                MALT10k = c(score_malt10k, non_score_malt10k, overall_score_malt10k),
                Golomb = c(score_zhang, non_score_zhang, overall_score_zhang))


df$algo = factor(rep(names(score_pbmc10k), 3), 
                 levels = c('raw', 'decontX', 'nb', 'poi', 'dsb', 'dsb_low', 'totalvi','scar'))
levels(df$algo) = c('raw', 'decontX','DecontPro', 'Neg.Bin.', 'dsb', 'dsb(high thre.)', 'totalVI', 'scAR')

df$group = factor(c(rep('Positive Score (Native Markers)', 8), 
                    rep('Negative Score (Other Cell Type Markers)',8),
                    rep('Overall Score',8)),
                    levels = c('Positive Score (Native Markers)',
                               'Negative Score (Other Cell Type Markers)',
                               'Overall Score'))

df.m = reshape2::melt(df)

## Paired percentage filling plot
df.m$bg = 1
F4B = ggplot(df.m, aes(variable, value, fill = group)) +
  geom_bar(stat = "identity", position="dodge") +
  geom_bar(aes(variable, bg),
           stat = "identity",
           position = "identity",
           fill = NA,
           color = "black") +
  scale_fill_jama(alpha = 0.7) +
  coord_cartesian(ylim = c(0, 1.05)) +
  # scale_fill_manual(values = c("#7CAE00", "#F8766D")) + 
  theme_bw() +
  labs(x = "", y = "Markers-Expressing Percentage", fill = "", tag = "B") + 
  theme(plot.title = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.margin = margin(-30,0,0,0),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9),
        strip.placement = "outside",
        plot.tag = element_text(face = "bold",
                                size = 15)) +
  geom_text(aes(y = bg, label = 100*round(value, 2)),
            position = position_dodge(0.9),
            vjust = -0.2,
            size = 2.8) +
  facet_wrap(~algo, strip.position = "top")



S17 = F4B + labs(tag = "")
ggsave("FigS17.pdf",
       S17,
       width = 11,
       height = 8.5,
       units = "in")




## ========== Density Plot

dataname = 'data_pbmc10khealthdonor'
dat = all_datasets[[dataname]]

adts = c('CD3', 'CD4', 'PD-1')

raw_df = data.frame(t(dat$raw[adts,]))
raw_df$source = 'Raw'
raw_df.m = reshape2::melt(raw_df)

# Poisson

poi_df = data.frame(t(dat$poi[adts,]))
poi_df$source = 'Decontaminated'
poi_df.m = reshape2::melt(poi_df)

poi_df.m = rbind(raw_df.m, poi_df.m)
poi_df.m$algo = 'DecontPro'



# dsb (need exponentiation)

dsb_df = 2^(data.frame(t(dat$dsb[adts,])))
dsb_df$source = 'Decontaminated'
dsb_df.m = reshape2::melt(dsb_df)

dsb_df.m = rbind(raw_df.m, dsb_df.m)
dsb_df.m$algo = 'dsb'


# totalVI

totalvi_df = data.frame(t(dat$totalvi[adts,]))
totalvi_df$source = 'Decontaminated'
totalvi_df.m = reshape2::melt(totalvi_df)

totalvi_df.m = rbind(raw_df.m, totalvi_df.m)
totalvi_df.m$algo = 'totalVI'


# scAR

scar_df = data.frame(t(dat$scar[adts,]))
scar_df$source = 'Decontaminated'
scar_df.m = reshape2::melt(scar_df)

scar_df.m = rbind(raw_df.m, scar_df.m)
scar_df.m$algo = 'scAR'


df.m = rbind(poi_df.m, dsb_df.m, totalvi_df.m, scar_df.m)
df.m$source = factor(df.m$source, levels = c('Raw', 'Decontaminated'))
df.m$algo = factor(df.m$algo, levels = c('DecontPro', 'dsb', 'totalVI', 'scAR'))


F4C = ggplot(df.m, aes(value, fill = source)) +
  geom_density(alpha = 0.5) +
  coord_cartesian(xlim = c(0, 10^5)) +
  scale_x_continuous(trans=scales::pseudo_log_trans(), breaks=c(10^(0:4)),
                     labels = function(x) format(x, scientific = FALSE)) +
  scale_fill_npg(alpha = 0.7) + 
  theme_bw() +
  theme(plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9),
        plot.margin = margin(0,2,0,0),
        legend.margin = margin(-30,0,0,0),
        legend.position = "bottom",
        plot.tag = element_text(face = "bold",
                                size = 15)) + 
  labs(x = "", y = "", fill ="PBMC10k", tag = 'C') +
  facet_grid(rows = vars(variable),
             cols = vars(algo),
             scales = "free")









## ============== Box plot
dataname = 'data_pbmc10khealthdonor'
dat = all_datasets[[dataname]]
cell_type = dat$adt_seurat$cellType


protein = c('PD-1')

df = data.frame(# CLR = dat$raw_norm[protein,], # counts(adt_filtered)[protein,],
                dsb = dat$dsb[protein,],
                scAR = dat$scar_denoised_norm[protein,], # scar_denoised[protein,],
                totalVI = dat$totalvi_denoised_norm[protein,], # totalvi_denoised[protein,],
                DecontPro = dat$poi_norm[protein,], # decontaminated_counts[protein,],
                cell_type = cell_type)





df.m = reshape2::melt(df)

df.m$variable = factor(df.m$variable, levels = c('DecontPro', 'dsb', 'totalVI', 'scAR'))


F4D = ggplot(df.m, aes(x = cell_type, y = value, fill = variable)) +
  geom_boxplot(lwd = 0.2, outlier.size = 0.5) +
  scale_y_continuous(trans=scales::pseudo_log_trans(), breaks=c(1,5,10^(1:4))) +
  scale_fill_npg(alpha = 0.7) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9),
        plot.margin = margin(0,0,0,0),
        legend.margin = margin(-30,0,0,0),
        legend.position = "bottom",
        plot.tag = element_text(face = "bold",
                                size = 15)) + 
  labs(x = "", y = "", fill = "PBMC10k", tag = "D") +
  ggtitle(paste0(protein, " after Decontamination and Scaling")) +
  facet_wrap(~variable,
             scales = "free")



## ====================
## ====== Figure 4
## ====================
# l = list()
# l[[1]] = F4A
# l[[2]] = F4B
# l[[3]] = F4C
# l[[4]] = F4D
# 
# 
# design = "
# AB
# CD
# "
# 
# F4 = wrap_plots(l) +
#    plot_layout(design = design) &
#   theme(legend.margin=margin(l = 0, r = 0))

F4 = (F4A + F4B) / (F4C + F4D)

ggsave("Figure4.pdf",
       F4,
       width = 11*1.3,
       height = 8.5*1.3,
       units = "in")




## ========== Ridge plot

# 
datadir_name = 'data_pbmc10khealthdonor'
adt_seurat = all_datasets$data_pbmc10khealthdonor$adt_seurat


datadir_name = 'data_pbmc5knextgem'
adt_seurat = all_datasets$data_pbmc5knextgem$adt_seurat


datadir_name = 'data_malt10k'
adt_seurat = all_datasets$data_malt10k$adt_seurat




# Import Ambient profile (see fig_modelApproximation.R)

df = read.csv(paste0("ambient_",datadir_name, '.csv'),
              row.names = 1)

ambient_profile = df$low
names(ambient_profile) = rownames(df)





## ridge plot 

proteins = c('CD3', 'CD4', "CD8a", "CD16", "CD19", "CD45RA")

counts = adt_seurat[['ADT']]@counts
counts_norm = sweep(counts, 2, colSums(counts), '/')

df_ridge <- data.frame(t(counts_norm[proteins,]))
df_ridge$cell_type = adt_seurat$cellType

df_ridge.m = reshape2::melt(df_ridge)



S9 = ggplot(df_ridge.m, aes(x = value, y = cell_type)) +
  geom_density_ridges(scale = 1, alpha = 0.5, fill = '#F8766D') +
  geom_vline(data = data.frame(variable = proteins,
                               ambient_profile_val = ambient_profile[proteins]),
             aes(xintercept = ambient_profile_val)) +
  theme_classic() +
  labs(x="", y="") +
  coord_cartesian(xlim = c(0,0.3)) +
  theme(legend.position="none",
        plot.tag = element_text(face = "bold",
                                size = 15)) +
  facet_wrap(~variable,
             strip.position = "bottom") +
  ggtitle('Library Size Normalized ADT Density')




ggsave(paste0("FigS9_",datadir_name,".pdf"), 
       S11,
       width = 11,
       height = 8.5,
       units = "in")





## === Compare Poisson and NB density plot (Supplementary)

## PBMC10k
## Density plot

density_plot = function(dat_list, title_text){
  
  proteins = rownames(dat_list$adt_filtered)
  p = list()
  
  for (i in 1:length(proteins)) {
    
    protein = proteins[i]
    df_poi <- data.frame(counts=assay(dat_list$adt_filtered,
                                      'counts')[protein,],
                         poi=dat_list$poi[protein,],
                         nb=dat_list$nb[protein,])
    df_poi.m = reshape2::melt(df_poi)
    
    # Poisson plot
    p1 = ggplot(df_poi.m, aes(value, fill = variable)) +
      geom_density(alpha = 0.5, adjust = 2) +
      scale_x_continuous(trans=scales::pseudo_log_trans(), breaks=c(1,5,10^(1:4))) +
      # scale_fill_discrete() +
      scale_fill_npg(alpha = 0.7, labels = c('Original', 'DecontPro', 'Neg.Bin.')) + 
      ggtitle(protein) +
      labs(x = "", fill =title_text) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.9),
            # legend.title = element_blank(),
            legend.margin=margin(t = -10))
    
    ylimit = layer_scales(p1)$y$get_limits()
    ylimit[2] = min(ylimit[2], 5)
    
    p[[i]] = p1 + coord_cartesian(ylim = ylimit)
  }
  
  return(p)
}

p = density_plot(all_datasets$data_pbmc10khealthdonor, "PBMC10k")

S18 = wrap_plots(p, ncol = 6, guides = "collect") +
  plot_layout()&theme(legend.position = "bottom", 
                      plot.margin = unit(c(3,3,2,1), "pt"))

ggsave("FigS18.pdf",
       S18,
       width = 11,
       height = 8.5,
       units = "in")



p = density_plot(all_datasets$data_pbmc5knextgem, "PBMC5k")

S19 = wrap_plots(p, ncol = 6, guides = "collect") +
  plot_layout()&theme(legend.position = "bottom", 
                      plot.margin = unit(c(3,3,2,1), "pt"))

ggsave("FigS19.pdf",
       S19,
       width = 11,
       height = 8.5,
       units = "in")



p = density_plot(all_datasets$data_malt10k, "MALT10k")

S20 = wrap_plots(p, ncol = 6, guides = "collect") +
  plot_layout()&theme(legend.position = "bottom", 
                      plot.margin = unit(c(3,3,2,1), "pt"))

ggsave("FigS20.pdf",
       S20,
       width = 11,
       height = 8.5,
       units = "in")



p = density_plot(all_datasets$young_aged, "Golomb")

S21 = wrap_plots(p, ncol = 6, guides = "collect") +
  plot_layout()&theme(legend.position = "bottom", 
                      plot.margin = unit(c(3,3,2,1), "pt"))

ggsave("FigS21.pdf",
       S21,
       width = 11,
       height = 8.5,
       units = "in")

