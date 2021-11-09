library(tidyverse)

rctd_results <- read.table("results_rctd_table.csv", sep = ',', header = TRUE, row.names = 1)

# Number of types above 30 beads (vline coordinate)
rctd_results %>% filter(spot_class %in% c('doublet_certain', 'doublet_uncertain', 'reject')) %>% {sum(count(., types)$n > 30)}

# Number of beads in types kept
rctd_results %>% filter(spot_class %in% c('doublet_certain', 'doublet_uncertain', 'reject')) %>% count(types) %>% filter(n > 30) %>% pull(n) %>% sum()

# Number of beads in types filtered out
rctd_results %>% filter(spot_class %in% c('doublet_certain', 'doublet_uncertain', 'reject')) %>% count(types) %>% filter(n <= 30) %>% pull(n) %>% sum()

rctd_hist <- rctd_results %>%
  filter(spot_class %in% c('doublet_certain', 'doublet_uncertain', 'reject')) %>%
  count(types) %>%
  arrange(desc(n)) %>%
  mutate(types = as_factor(types),
         included = n > 30) %>%
  ggplot() +
  geom_bar(aes(types, n), size = 0, stat = 'identity', col = 'gray30', fill = 'gray30') + 
  # geom_hline(yintercept = 30) + 
  geom_vline(xintercept = 84, size = 1, linetype = 'longdash', alpha = 0.5) +
  theme_minimal() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid = element_blank()) +
  labs(x = 'Doublet type', y = 'Number of beads') + 
  annotate('text', x = 50, y = 10300, label = glue::glue('16/26: 10797 beads'), color = 'gray20') +
  geom_curve(
    aes(x = 37, y = 9950, xend = 3, yend = 9250),
    arrow = arrow(length = unit(0.07, "inch"), type = 'closed'), size = 0.5,
    color = "gray20", curvature = -0.25
  ) + 
  annotate('text', x = 110, y = 1250, label = '3/31: 30 beads', color = 'gray20') +
  geom_curve(
    aes(x = 102, y = 900, xend = 85, yend = 100),
    arrow = arrow(length = unit(0.07, "inch"), type = 'closed'), size = 0.5,
    color = "gray20", curvature = -0.05
  )
ggsave('rctd_type_hist.png', rctd_hist, width = 7.5, height = 3, units = 'in')

all_markers <- read_csv('all_markers.csv')

# marker_pct <- all_markers %>%
#   mutate(diff = pct.2 - pct.1) %>%
#   pivot_longer(c(pct.1, pct.2), names_to = 'pct.type', values_to = 'pct') %>%
#   mutate(ID = tidytext::reorder_within(ID, diff, types)) %>%
#   ggplot() +
#   geom_col(aes(x = ID, y = pct, fill = pct.type), position = 'dodge') +
#   facet_wrap(~types, scales = 'free_x', nrow = 6) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90),
#         strip.background = element_rect(fill = 'white'),
#         legend.position = c(0.865, 0.08),
#         legend.title.align = 0.5) +
#   ylab('Detection fraction') +
#   scale_fill_discrete(name = 'Type', labels = c('Type 1', 'Type 2')) +
#   tidytext::scale_x_reordered()

marker_pct <- all_markers %>%
  mutate(diff = pct.2 - pct.1) %>%
  pivot_longer(c(pct.1, pct.2), names_to = 'pct.type', values_to = 'pct') %>%
  mutate(ID = tidytext::reorder_within(ID, diff, types)) %>%
  ggplot() + 
  geom_col(aes(x = ID, y = pct, fill = pct.type), position = 'dodge') + 
  facet_wrap(~types, scales = 'free_y', nrow = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25),
        strip.background = element_rect(fill = 'white'),
        legend.position = c(0.96, 0.15),
        legend.title.align = 0.5) +
  labs(x = NULL, y = NULL) +
  scale_fill_discrete(name = 'Type', labels = c('First type', 'Second type')) + 
  coord_flip() + 
  tidytext::scale_x_reordered()
ggsave('marker_pct_barplot.png', marker_pct, width = 11.3, height = 10, units = 'in')
