library(RCTD)
library(Matrix)
library(dplyr)

refdir <- ""
spatialdir <- ""

# Construct a reference based on single cell profiles
ref_counts <- readRDS(file.path(refdir, 'scRNASeq_common_genes_counts.RDS'))
ref_meta <- read.csv(file.path(refdir, 'meta_data_common_genes.csv'), header = TRUE)
ref_types <- ref_meta$cluster %>% set_names(ref_meta$barcode) %>% as.factor()
ref_numi <- ref_meta$nUMI %>% set_names(ref_meta$barcode)
ref <- Reference(counts = ref_counts, cell_types = ref_types, nUMI = ref_numi)

# Construct a puck based on the spatial transcriptomics reads
st_counts <- read.csv(file.path(spatialdir, 'MappedDGEForR_common_genes.csv'), row.names = 1)
st_coords <- read.csv(file.path(spatialdir, 'BeadLocationsForR_common_genes.csv'), row.names = 1)
st_numi <- colSums(st_counts)
puck <- SpatialRNA(st_coords, st_counts, st_numi)

# Run RCTD
myRCTD <- create.RCTD(puck, ref, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

# Save full RDS object
saveRDS(myRCTD, 'velocity_files/results_rctd_fulldata.rds')
print('Saved results_rctd_fulldata.rds.')

# Save results table
classes <- myRCTD@results$results_df
classes$barcode <- rownames(classes)

proportions <- myRCTD@results$weights_doublet %>% 
    as.matrix() %>%
    as.data.frame() %>%
    rename(first_type_prop = first_type, second_type_prop = second_type)
proportions$barcode <- rownames(proportions)

doublet_classes <- c('reject', 'doublet_certain', 'doublet_uncertain')

results <- full_join(classes, proportions) %>%
    rowwise() %>% 
    mutate(type1 = if_else(spot_class %in% doublet_classes,
                            min(as.numeric(first_type), as.numeric(second_type)),
                            as.numeric(first_type)),
            type2 = if_else(spot_class %in% doublet_classes,
                            max(as.numeric(first_type), as.numeric(second_type)),
                            as.numeric(first_type)),
            types = if_else(type1 != type2,
                            paste(type1, type2, sep = '/'),
                            as.character(type1)),
            prop1 = if_else(spot_class %in% doublet_classes & as.numeric(first_type) > as.numeric(second_type),
                            second_type_prop,
                            first_type_prop),
            prop2 = if_else(spot_class %in% doublet_classes & as.numeric(first_type) > as.numeric(second_type),
                            first_type_prop,
                            second_type_prop)
            ) %>% 
    ungroup() %>%
    select(barcode, spot_class, type1, type2, types, prop1, prop2, first_type, second_type, first_type_prop, second_type_prop, min_score, singlet_score)

write.csv(results, 'velocity_files/results_rctd_table.csv')
print('Saved results_rctd_table.csv. Script finished!')
