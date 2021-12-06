# Small script to compare the detected markers to the pre-found markers

library(dplyr)
library(readr)
library(tidyr)

cell_markers <- read_csv('XXX')
found_markers <- read_csv('XXX')

cell_markers_long <- pivot_longer(cell_markers, -c(Cluster, Celltype, Abbrev), names_to = 'n_marker', values_to = 'ID')

# Print marker overlaps
print('Shared markers:')
intersect(cell_markers_long$ID, found_markers$ID)
print('Unique to cell markers:')
setdiff(cell_markers_long$ID, intersect(cell_markers_long$ID, found_markers$ID))
print('Unique to detected markers:')
setdiff(found_markers$ID, intersect(cell_markers_long$ID, found_markers$ID))

# Get shared marker information
mrk_info <- inner_join(cell_markers_long, found_markers)
print(mrk_info)
print(count(mrk_info, Cluster, Celltype))


joined_markers <- bind_rows(
    cell_type_markers = select(mutate(cell_markers_long, type = as.character(Cluster)), type, ID), 
    additional_markers = select(found_markers, type = types, ID), 
    .id = 'origin') %>% 
    filter(ID != '-')
write_csv(joined_markers, 'velocity_files/joint_cell_type_markers.csv')

