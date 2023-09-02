params <- commandArgs(T)
stanfun <- params[1]
datadir_name <- params[2]


#############################################################################


library(celda)
library(DropletUtils)
library(ggplot2)
library(patchwork)
library(ggridges)
library(Seurat)
library(rstan)
rstan_options(auto_write = TRUE)

set.seed(12345)


datadir = '/rprojectnb2/camplab/home/yin/decontX_ADT/zhangetal/young-aged/'




# Import filtered
sce_filtered = read10xCounts(paste0(datadir,'/filtered_feature_bc_matrix'))
sce_raw = read10xCounts(paste0(datadir,'/raw_feature_bc_matrix'))


# Set rowname and colname
rownames(sce_filtered) = rowData(sce_filtered)$Symbol
colnames(sce_filtered) = colData(sce_filtered)$Barcode

rownames(sce_raw) = rowData(sce_raw)$Symbol
colnames(sce_raw) = colData(sce_raw)$Barcode


gex_filtered = sce_filtered[rowData(sce_filtered)$Type == 'Gene Expression',]
gex_raw = sce_raw[rowData(sce_raw)$Type == 'Gene Expression',]

adt_filtered = sce_filtered[rowData(sce_filtered)$Type == 'Antibody Capture',]
adt_raw = sce_raw[rowData(sce_raw)$Type == 'Antibody Capture',]


# rownames(adt_raw) = rowData(adt_raw)$ID # Remove the annoying "_TotalSeqB" suffix



# HTO filter
hto_f = grepl('HTO', rownames(adt_raw))

hto_filtered = adt_filtered[hto_f, ]
hto_raw = adt_raw[hto_f, ]

adt_filtered = adt_filtered[!hto_f, ]
adt_raw = adt_raw[!hto_f, ]


# Quick diagnostic plot
barplot(apply(counts(hto_filtered), 1, max))



# Filter out low expression ADT
adt_f = rowSums2(counts(adt_filtered) > 0) > 50
adt_filtered = adt_filtered[adt_f, ]
adt_raw = adt_raw[adt_f, ]



rownames(adt_filtered) = substring(rownames(adt_filtered), 6)
rownames(adt_raw) = substring(rownames(adt_raw), 6)



## ==== QC

# Set filter
mtgene = c(grep(pattern = "^mt-", rownames(gex_raw), value = TRUE),
           grep(pattern = "^MT-", rownames(gex_raw), value = TRUE))

md = data.frame(
  rna.size = colSums2(counts(gex_raw)), 
  prot.size = colSums2(counts(adt_raw)), 
  n.gene = colSums2(counts(gex_raw) > 0), 
  mt.prop = colSums2(counts(gex_raw[mtgene, ])) / colSums2(counts(gex_raw))
)

rownames(md) = colnames(sce_raw)


md$drop.class = ifelse(rownames(md) %in% colnames(sce_filtered), 'cell', 'background')
cellmd = md[md$drop.class == 'cell', ]



# Remove 0 count drop
md = md[md$rna.size > 0 & md$prot.size > 0, ]

# QC background
background_drops = rownames(
  md[md$rna.size < 100 &
       md$prot.size <25, ]
) 


adt_empty = adt_raw[,colnames(adt_raw) %in% background_drops]



# QC cell

qc_cells = rownames(
  cellmd[cellmd$rna.size > 0 &
           cellmd$prot.size > 0 &
           cellmd$rna.size > quantile(cellmd$rna.size, 0.01) &
           cellmd$rna.size < quantile(cellmd$rna.size, 0.99) &
           cellmd$prot.size > quantile(cellmd$prot.size, 0.01) &
           cellmd$prot.size < quantile(cellmd$prot.size, 0.99) &
           cellmd$mt.prop < 0.14,]
)


adt_filtered = adt_raw[, colnames(adt_raw) %in% qc_cells]
gex_filtered = gex_raw[, colnames(gex_raw) %in% qc_cells]





# Find isotypes
isotype = rownames(adt_filtered)[grepl('IgG', rownames(adt_filtered))]





# save original adt_filtered
adt_w_isotype = adt_filtered


adt_filtered = adt_filtered[!(rownames(adt_filtered) %in% isotype),]
adt_empty = adt_empty[!(rownames(adt_empty) %in% isotype),]




# Seurat to get clusters
adt_seurat = Seurat::CreateSeuratObject(counts(adt_filtered),
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
                           resolution = 0.3,
                           verbose = FALSE)

adt_seurat <- RunUMAP(adt_seurat,
                      dims = 1:10,
                      assay = "ADT",
                      reduction = "pca_adt",
                      reduction.name = "adtUMAP",
                      verbose = FALSE)

adt_seurat[["adtClusterID"]] <- Idents(adt_seurat)





cell_type = adt_seurat[["adtClusterID"]]
cell_type = as.integer(cell_type$adtClusterID)


# DecontX
adt_filtered = decontX(adt_filtered, z = cell_type)



# Contamination rate p

counts = as.matrix(counts(adt_filtered))

p = rowSums2(counts(adt_empty))
p = p/sum(p)

# P from cells
p_cell = rowSums2(counts(adt_filtered))
p_cell = p_cell/sum(p_cell)



## ======= rstan

dat = list(N = nrow(counts),
           M = ncol(counts),
           K = length(unique(cell_type)),
           cell_type = cell_type,
           counts = counts,
           OC = colSums(counts),
           run_estimation = 1,
           p = p_cell,
           delta_sd = 2e-5,
           background_sd = 2e-6)



m = stan_model(file = paste0(stanfun,'.stan'))

out = vb(object = m,
         init = list(delta = matrix(rep(1e-4, ncol(counts)*nrow(counts)),
                                    nrow = ncol(counts),
                                    ncol = nrow(counts)),
                     background = matrix(rep(1e-2, ncol(counts)*nrow(counts)),
                                         nrow = nrow(counts),
                                         ncol = ncol(counts))),
         data = dat,
         seed = 12345,
         iter = 50000)

save(dat, 
     out, 
     adt_filtered,
     gex_filtered,
     adt_empty,
     adt_seurat,
     adt_w_isotype,
     p,
     p_cell,
     file = paste0(stanfun,'_',datadir_name,'_Robj.Rdata'))

rm(list = ls())


