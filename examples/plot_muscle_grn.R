library(monocle3)
library(msigdbr)
library(igraph)
library(tidygraph)
library(network)
library(ggnetwork)
library(tibble)
library(googledrive)


drive_download("https://drive.google.com/file/d/1GD7fkEwVV9pclp1rYG5ygpFw7WMp321o/")

g <- readRDS("muscle_ggnetwork.rds")

colors = c( "#735CDD", "#1B9AAA","#EF476F", "#101D42")
names(colors) = c("mesodermal progenitor cells (contains PSM)",
                  "myoblast/fast-committed myocyte",
                  "fast-committed myocyte",
                  "mature fast muscle")

# lower cap sizing, otherwise doesn't show text
g = g %>% mutate(regulator_score = ifelse(regulator_score < 0.55, 0.55, regulator_score))


p = ggplot(mapping = aes(x, y, xend = xend, yend = yend)) +
  # draw undirected edges
  # draw activator edges
  geom_nodes(data = g, size=0.5) +
  geom_edges(data = g %>% filter(edge_type == "activator"),
             aes(size = 0.35*scaled_weight),
             arrow = arrow(length = unit(1.5, "pt"), type="closed"),
             curvature = 0.1) +
  # draw repressor edges
  geom_edges(data = g %>% filter(edge_type == "repressor"),
             aes(size = 0.35*scaled_weight),
             arrow = arrow(angle = 90, length = unit(.075, "cm")),
             curvature = 0) +
  geom_nodelabel_repel(data = g ,
                       aes(fill = cell_type_broad,
                           label = gene_short_name,
                           size = 1.25*abs(regulator_score)),
                       color="white",
                       label.size = 0.2,
                       #box.padding = 0.1,
                       box.padding = unit(0.1, "lines"),
                       label.padding = 0.15,
                       force = 0.1) +
  scale_fill_manual(values = colors) +
  scale_size_identity() +
  monocle3:::monocle_theme_opts() +
  theme_blank() +
  theme(legend.position = "none")

ggsave(p, filename = "muscle_grn.png", height = 4, width = 3)
