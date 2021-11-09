# Load libraries and data --------------------------------------

library(RCTD)
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)
library(purrr)
library(readr)
library(Seurat)
library(glue)
library(tidyr)
library(ggtext)

RCTD_path <- ""
slideseq_path <- ""
celltypes_path <- ""
locations_path <- "" 

if (!dir.exists(file.path('heatmaps'))) {
  dir.create(file.path('heatmaps'))
}
if (!dir.exists(file.path('velocity_files'))) {
  dir.create(file.path('velocity_files'))
}

add_meta <- read.table(RCTD_path, sep = ',', header = TRUE, row.names = 1)

# Gather data in Seurat object --------------------------------------
seur <- read_table(slideseq_path)
seur <- column_to_rownames(seur, 'GENE')

# Create the object and add spatial data
seur <- CreateSeuratObject(seur, assay = 'Spatial')
seur@images$image <- ReadSlideSeq(locations_path)
seur@meta.data <- rownames_to_column(seur@meta.data, 'barcode') %>%
  left_join(add_meta, by = 'barcode') %>%
  column_to_rownames('barcode')

# Normalise and scale counts
seur <- seur %>% 
  subset(subset = nCount_Spatial > 200) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)
seur <- ScaleData(seur, assay = 'Spatial')
Idents(seur) <- "types"

# Show and visualise distributions --------------------------------------
# Update to take combinations filtering into account

# Print number of singlets/doublets
print("Spot classes of QC'd spots:")
seur@meta.data %>% 
  count(spot_class) %>% 
  print()

# Print type distributions
print("Top type combinations of QC'd spots:")
seur@meta.data %>% 
  count(types) %>% 
  filter(n > 30) %>% 
  print()

print("Total number of spots after QC:")
print(nrow(seur@meta.data))

  

# # Plot spatial distributions
# # Show multiple stats at once? 
# BetterSpatialPlot <- function(seurat_obj, feat_name) {
#   coords <- seurat_obj@images$image@coordinates
#   feat <- seurat_obj@meta.data %>% 
#     rownames_to_column('barcode') %>% 
#     rename(feature = !!sym(feat_name))
#   plot_data <- full_join(x = coords, y = feat, by = c('cells' = 'barcode'))
#   plot_data <- plot_data %>% mutate(feature = fct_lump_n(feature, 7))
#   # plot_data <- plot_data %>% filter(feature %in% c('16', '16/26', '26'))
#   ggplot(plot_data) + 
#     geom_point(aes(x, y, col = feature), size = 0.5) + 
#     guides(col = guide_legend(title = feat_name)) +
#     theme_minimal()
# }

# BetterSpatialPlot(seur, 'types')
# BetterSpatialPlot(seur, 'spot_class')

# Analyse top genes --------------------------------------

over_3 <- seur@meta.data %>% count(types) %>% filter(n > 3) %>% pull(types)

combinations <- seur@meta.data %>% 
  # Get counts of the doublets
  filter(spot_class %in% c('reject', 'doublet_certain', 'doublet_uncertain')) %>%
  count(type1, type2) %>%
  # Keep only doublets with over 30 beads
  filter(n > 30) %>%
  mutate(types = paste0(type1, '/', type2)) %>%
  # Keep only doublets whose matching singlets at least 3 beads
  filter(type1 %in% over_3, type2 %in% over_3) %>%
  arrange(desc(n)) %>%
  select(-n)

combinations_list <- setNames(
  pmap(combinations, ~c(type1 = ..1, types = ..3, type2 = ..2)),
  combinations$types
) 


# Map over all comparisons to find markers and build heatmap
all_markers <- map(combinations_list, function(combs) {
  # Find markers 
  print(glue("Finding markers for {combs['type1']}/{combs['type2']}"))
  mrks <- FindMarkers(seur, ident.1 = combs['type1'], ident.2 = combs['type2'], min.diff.pct = 0.67) %>% 
    filter(pct.1 < 0.1 | pct.2 < 0.1)
  print(glue('Number of markers: {nrow(mrks)}'))
  seur_combs <- seur %>% subset(types %in% combs)

  if (nrow(mrks) > 1) {    
    zscores <- map_dfr(combs, 
      function(types_value) {
        data <- seur_combs %>%
          subset(types == types_value) %>% 
          GetAssayData('scale.data') %>%
          .[rownames(mrks),]
        
        ordering <- hclust(dist(t(data)), method = "average")$order
        # ordering <- get_order(seriate(dist(scale(t(data))), 'HC')) #requires library(seriation)
      
        data %>% 
          as.data.frame() %>% 
          rownames_to_column('marker') %>% 
          pivot_longer(cols = -marker, names_to = 'barcode', values_to = 'zscore') %>%
          mutate(
            facet_text = as.factor(glue("{types_value}<br><span style='font-size:9pt'>{length(unique(barcode))} beads, {length(unique(marker))} markers</span>")),
            marker = factor(marker, levels = rownames(arrange(mrks, pct.1 - pct.2))),
            barcode = factor(barcode, levels = barcode[ordering]))
      })
    
    
    # v4
    heat <- zscores %>% 
      ggplot(aes(x = barcode, y = marker, fill = zscore)) + 
        geom_tile() +
        theme(axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_blank(),
              strip.text.x = ggtext::element_markdown(size = 12),
              panel.spacing = unit(1, "lines")
              ) +
        scale_fill_distiller(name = "z-scores", 
                            # limits = c(-1, 4),
                            palette = 'RdBu') +
        facet_wrap(~facet_text, scales = 'free_x')
      # ggsave('test_heatmap2.png', heat)

    ggsave(glue("heatmaps/ggheatmap2_{combs['type1']}_{combs['type2']}.png"), heat)
    print(glue("Finished plot ggheatmap2_{combs['type1']}_{combs['type2']}.png"))

    # Return markers as output to list
    return(arrange(mrks, desc(pct.1 - pct.2)))
  } else {
    return(data.frame(p_val = NA, avg_log2FC = NA, pct.1 = NA, pct.2 = NA, p_val_adj = NA)) }
})

# Create output objects
markers_output <- all_markers %>% 
  map(~rownames_to_column(.x, 'ID')) %>%
  bind_rows(.id = 'types') %>%
  filter(!is.na(p_val))

# Write output to files
write_csv(markers_output, 'velocity_files/all_markers.csv')

print(glue("all_markers.csv: {nrow(markers_output)} type/marker combinations"))
