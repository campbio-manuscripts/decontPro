library(DropletUtils)
library(dsb)
library(ggplot2)
library(patchwork)
library(scattermore)



datadir_name = 'data_pbmc10khealthdonor'
# datadir_name = 'data_pbmc5knextgem'
# datadir_name = 'data_malt10k'
datadir = paste0('/rprojectnb2/camplab/home/yin/poisson/', datadir_name)


# datadir_name = 'young_aged'
# datadir = '/rprojectnb2/camplab/home/yin/decontX_ADT/zhangetal/young-aged'


# Import
sce_filtered = read10xCounts(paste0(datadir,'/filtered_feature_bc_matrix'))
colnames(sce_filtered) = colData(sce_filtered)$Barcode

sce_raw = read10xCounts(paste0(datadir,'/raw_feature_bc_matrix'))
colnames(sce_raw) = colData(sce_raw)$Barcode


gex_filtered = sce_filtered[rowData(sce_filtered)$Type == 'Gene Expression',]
gex_raw = sce_raw[rowData(sce_raw)$Type == 'Gene Expression',]

# Assign gene name
rownames(gex_filtered) = rowData(gex_filtered)$Symbol
rownames(gex_raw) = rowData(gex_raw)$Symbol


adt_filtered = sce_filtered[rowData(sce_filtered)$Type == 'Antibody Capture',]
adt_raw = sce_raw[rowData(sce_raw)$Type == 'Antibody Capture',]


rownames(adt_raw) = rowData(adt_raw)$ID # Remove the annoying "_TotalSeqB" suffix



# HTO filter
hto_f = grepl('HTO', rownames(adt_raw))
adt_filtered = adt_filtered[!hto_f, ]
adt_raw = adt_raw[!hto_f, ]

hto_f = grepl('^ADT', rownames(adt_raw))
adt_filtered = adt_filtered[!hto_f, ]
adt_raw = adt_raw[!hto_f, ]


# Filter out low expression ADT
adt_f = rowSums2(counts(adt_filtered) > 0) > 50
adt_filtered = adt_filtered[adt_f, ]
adt_raw = adt_raw[adt_f, ]







# Set filter
mtgene = c(grep(pattern = "^mt-", rownames(gex_raw), value = TRUE),
           grep(pattern = "^MT-", rownames(gex_raw), value = TRUE))

md = data.frame(
  rna.size = colSums2(counts(gex_raw)), 
  prot.size = colSums2(counts(adt_raw)), 
  n.gene = colSums2(counts(gex_raw) > 0), 
  mt.prop = colSums2(counts(gex_raw[mtgene, ])) / colSums2(counts(gex_raw))
)

md$rna.size.log = log10(md$rna.size)
md$prot.size.log = log10(md$prot.size)

rownames(md) = colnames(sce_raw)


md$drop.class = ifelse(rownames(md) %in% colnames(sce_filtered), 'cell', 'background')
cellmd = md[md$drop.class == 'cell', ]


# Remove 0 count drop
md = md[md$rna.size > 0 & md$prot.size > 0, ]


# Check by plotting
scattermoreplot(md$rna.size.log,
                md$prot.size.log,
                # col = as.integer(md$drop.class) + 1,
                cex=0.5,
                xlab = 'RNA Lib Size',
                ylab = 'ADT Lib Size')
# abline(h = log10(500))


# QC background
background_drops = rownames(
  md[ md$prot.size.log > 1.5 & 
        md$prot.size.log < 3 & 
        md$rna.size.log < 2.5, ]
) 


background.adt.mtx = as.matrix(counts(adt_raw[,colnames(adt_raw) %in% background_drops]))



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


cell.adt.mtx = as.matrix(counts(adt_raw[, colnames(adt_raw) %in% qc_cells]))




isotype.controls = rownames(cell.adt.mtx)[grepl('IgG', rownames(cell.adt.mtx))]


# Denoise with dsb
if (identical(isotype.controls, character(0))) {
  cell.dsb.norm = DSBNormalizeProtein(
    cell_protein_matrix = cell.adt.mtx, 
    empty_drop_matrix = background.adt.mtx, 
    denoise.counts = TRUE, 
    use.isotype.control = FALSE
  )
} else {
  cell.dsb.norm = DSBNormalizeProtein(
    cell_protein_matrix = cell.adt.mtx, 
    empty_drop_matrix = background.adt.mtx, 
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, 
    isotype.control.name.vec = isotype.controls
    # use.isotype.control = FALSE
  )
}



# Save results
save(cell.adt.mtx,
     cell.dsb.norm,
     file = paste0(datadir_name,'_dsb.Rdata'))
