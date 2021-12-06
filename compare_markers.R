myPaths <- c("/home/cblaauw/library", .libPaths())
.libPaths(myPaths)

library(dplyr)
library(readr)
library(tidyr)

cell_markers <- read_csv('project_data/scRNA-seq/embryo_cell_type_markers.csv')
all_markers_050 <- read_csv('velocity_files/all_markers_050.csv')

cell_markers_long <- pivot_longer(cell_markers, -c(Cluster, Celltype, Abbrev), names_to = 'n_marker', values_to = 'ID')

# Print marker overlaps
print('Shared markers:')
intersect(cell_markers_long$ID, all_markers_050$ID)
print('Unique to cell markers:')
setdiff(cell_markers_long$ID, intersect(cell_markers_long$ID, all_markers_050$ID))
print('Unique to detected markers:')
setdiff(all_markers_050$ID, intersect(cell_markers_long$ID, all_markers_050$ID))

# Get shared marker information
mrk_info <- inner_join(cell_markers_long, all_markers_050)
print(mrk_info)
print(count(mrk_info, Cluster, Celltype))


joined_markers <- bind_rows(
    cell_type_markers = select(mutate(cell_markers_long, type = as.character(Cluster)), type, ID), 
    additional_markers = select(all_markers_050, type = types, ID), 
    .id = 'origin') %>% 
    filter(ID != '-')
write_csv(joined_markers, 'velocity_files/joint_cell_type_markers.csv')

