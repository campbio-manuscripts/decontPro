library(rstan)
library(bayesplot)
library(celda)
library(ggplot2)
library(patchwork)
library(Seurat)
library(scattermore)
library(ggridges)



# file_name = 'shrinkage_data_pbmc5knextgem'
# file_name = 'shrinkage_data_malt10k'
# file_name = 'shrinkage_data_pbmc10khealthdonor'


load(paste0(file_name,'_Robj.Rdata'))

# Extract optimized param values
val = out@sim$est


r_est = val$r


# If delta a matrix
delta_est = t(val$delta)
delta_mean_est = val$delta_mean

background_est = val$background
background_mean_est = val$background_mean



# =========== 

# 
p_est = dat$p


# If p fixed
unscaled_rates = matrix(p_est,
                        nrow = nrow(dat$counts),
                        ncol = ncol(dat$counts))

# If delta a matrix
scaling_factor = matrix(dat$OC,
                        nrow = nrow(dat$counts),
                        ncol = length(dat$OC),
                        byrow = T)

empty_rate_est = scaling_factor * unscaled_rates * delta_est * (1 - background_est)



unscaled_rates = t(r_est[dat$cell_type,])

# If delta a matrix
scaling_factor = matrix(dat$OC,
                        nrow = nrow(dat$counts),
                        ncol = length(dat$OC),
                        byrow = T)

cell_rate_est = scaling_factor * unscaled_rates * (1-delta_est) * (1 - background_est)


background_est = val$background

scaling_factor = matrix(dat$OC,
                        nrow = nrow(dat$counts),
                        ncol = length(dat$OC),
                        byrow = T)

background_rate_est = background_est * scaling_factor

rownames(background_est) = rownames(adt_filtered)
rownames(background_rate_est) = rownames(adt_filtered)


# Decontaminated counts
counts = dat$counts

#
decontaminated_counts = counts * cell_rate_est/(empty_rate_est + cell_rate_est + background_rate_est)
background_counts = counts * background_rate_est/(empty_rate_est + cell_rate_est + background_rate_est)

mean_param = empty_rate_est + cell_rate_est + background_rate_est




save(decontaminated_counts, 
     adt_filtered, 
     gex_filtered,
     adt_seurat, 
     file = paste0(file_name,'_poisson.Rdata'))
