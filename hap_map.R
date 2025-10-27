#setwd()

if (!require(ape)) install.packages("ape")
if (!require(devtools)) install.packages("devtools")
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggnewscale)) install.packages("ggnewscale")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggrepel)) install.packages("ggrepel")
if (!require(pegas)) install.packages("pegas")
if (!require(tidyr)) install.packages("tidyr")
if (!require(scatterpie)) devtools::install_github("YuLab-SMU/scatterpie")

library(ape)
library(devtools)
library(dplyr)
library(ggnewscale)
library(ggplot2)
library(ggrepel)
library(pegas)
library(scatterpie)
library(tidyr)

#############
# FUNCTIONS #
#############

prepare_haplo_edges <- function(net, hap_coords) {
  edge_mat <- unclass(net)                 
  labels <- attr(net, "labels")            
  
  edges <- data.frame(
    from   = labels[edge_mat[, 1]],
    to     = labels[edge_mat[, 2]],
    weight = edge_mat[, 3],   
    prob   = edge_mat[, 4],   
    stringsAsFactors = FALSE
  )
  
  edges <- edges %>%
    left_join(hap_coords[, c("Haplotype","Longitude","Latitude")],
              by = c("from" = "Haplotype")) %>%
    rename(x = Longitude, y = Latitude) %>%
    left_join(hap_coords[, c("Haplotype","Longitude","Latitude","Elevation")],
              by = c("to" = "Haplotype")) %>%
    rename(xend = Longitude, yend = Latitude, elev = Elevation)
  
  return(edges)
}

prep_hap_nodes <- function(hap_coords) {
  #make pie cols for populations
  pop_pie_cols <- setdiff(colnames(hap_coords), c("Haplotype", "Latitude", "Longitude", "Elevation"))
  #make pie cols for haplotypes
  hap_pie_cols <- unique(hap_coords$Haplotype)
  list(
    nodes = hap_coords,
    pop_pie_cols = pop_pie_cols,
    hap_pie_cols = hap_pie_cols
  )
}

#############
# INPUT DATA#
#############

groups <- list(
  am = 1:5,
  cl = 6:8,
  ec = 9:15,
  em = 16:21,
  et = 22:26,
  mt = 27:30,
  ne = 31:39,
  ty = 40:43
)
seq <- read.dna("haplo_co1.fasta", format="fasta")
pop.assign <- read.csv("haplo_pop.csv", header = TRUE, stringsAsFactors = FALSE) 

################################
# From Latest_sampling_sites.R #
################################
source("samp_site.R")

####################
# DATA PREPARATION #
####################

seq_list <- lapply(groups, function(idx) seq[idx, ])
pop_list <- lapply(groups, function(idx) as.character(pop.assign[idx, 2]))
hap_list <- lapply(seq_list, function(hp) pegas::haplotype(hp, strict = T))
hname_list <- lapply(hap_list, function(hn) paste("H", 1 : nrow(hn), sep = ""))
hap_list <- lapply(hap_list, function(h) {
  rownames(h) <- paste0("H", seq_len(nrow(h)))
  h
})
net_list <- lapply(hap_list, function(n) haploNet(n))
hap_table_list <- mapply(function(h, pops, seqs) {
  ind2hap <- attr(h, "index")
  hap_df <- do.call(rbind, lapply(seq_along(ind2hap), function(i) {
    data.frame(
      Taxa = rownames(seqs)[ind2hap[[i]]],
      Haplotype = paste0("H", i),
      Pop = pops[ind2hap[[i]]],
      stringsAsFactors = FALSE
    )
  }))
  pop_levels <- unique(pops)
  pop_mat <- sapply(pop_levels, function(p) as.integer(hap_df$Pop == p))
  colnames(pop_mat) <- pop_levels
  df <- cbind(hap_df[, c("Taxa", "Haplotype")], pop_mat)
  rownames(df) <- df$Taxa
  df$Taxa <- NULL
  df <- df[, c("Haplotype", pop_levels)]
  df
}, hap_list, pop_list, seq_list, SIMPLIFY = FALSE)

pop_list_extracted <- lapply(hap_table_list, function(hap.table) {
  colnames(hap.table)[-1]
})

pop.col <- c("#E69F00", "#56B4E9", "#F0E442", 
             "#009E73", "#CC79A7", "#0072B2", 
             "#920000", "#BDBDBD", "#AE98FF",
             "#EFEFEF", "#885D8A", "#D55E00"
             )

pop.level <- unique(pop.assign$Population)
hap.type <- unique(unlist(lapply(hap_table_list, function(hap) hap$Haplotype)))
pop_assign_col <- setNames(pop.col[seq_along(pop.level)], pop.level)
hap_assign_col <- setNames(pop.col[seq_along(hap.type)], hap.type)
hap_coords_list <- lapply(names(hap_table_list), function(grp) {
  hap_tab <- hap_table_list[[grp]]
  hap_names <- hap_tab$Haplotype
  taxa_names <- rownames(hap_tab)
  pop_names <- setdiff(colnames(hap_tab), c("Haplotype", "Taxa"))
  hap_mat <- hap_tab[, pop_names, drop = FALSE]
  coords <- t(sapply(taxa_names, function(taxon) {
    idx <- match(taxon, pop.assign$Taxa)
    return(c(pop.assign$Latitude[idx], 
               pop.assign$Longitude[idx], 
               pop.assign$Elevation[idx]))
  }))
  colnames(coords) <- c("Latitude","Longitude","Elevation")
  df <- data.frame(
    Haplotype = hap_names,
    coords,
    hap_mat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(df) <- taxa_names
  df
})
names(hap_coords_list) <- names(hap_table_list)
grp <- "em"
net <- net_list[[grp]]
hap_coords <- hap_coords_list[[grp]]
edges_df <- prepare_haplo_edges(net, hap_coords)
nodes_info <- prep_hap_nodes(hap_coords)
hap_scaled <- hap_coords %>%
  mutate(
    fake_long = vert_bound + (Elevation / max_elev)
  )
pop_cols <- setdiff(colnames(hap_scaled), c("Haplotype", "Latitude", "Longitude", "Elevation", "fake_long"))
hap_points <- hap_scaled %>%
  dplyr::mutate(
    Population = apply(
      dplyr::select(., all_of(pop_cols)),
      1,
      function(x) names(x)[which.max(x)]
    )
  )
hap_points_reshaped <- hap_points %>%
  dplyr::mutate(value = 1, id = dplyr::row_number()) %>%
  tidyr::pivot_wider(
    names_from = Haplotype,
    values_from = value,
    values_fill = list(value = 0)
  ) %>%
  dplyr::relocate(id, .before = 1)
midpoints <- edges_df %>%
  rowwise() %>%
  mutate(mx = mean(c(x, xend), na.rm = TRUE),
         my = mean(c(y, yend), na.rm = TRUE)
         )

########
# Plot #
########

main_haplo <- ggplot() +
  # === Elevation base ===
  geom_spatraster(data = elv_main_masked) + #, maxcell = 3.8e7) +
  scale_fill_gradientn(
    colours = alpha(alt_colors(100), 0.7), 
    na.value = NA, 
    guide = guide_colorbar(ticks.colour = "#000000", order = 1)
  ) +
  # Taiwan main island border
  geom_sf(data = tw_isle, color = "#000000", fill = NA) +
  ggnewscale::new_scale_fill() +
  # Elevation profile plots
  geom_scatterpie(data = hap_points_reshaped, 
                  aes(x = fake_long, y = Latitude, r = 0.04), 
                  cols = nodes_info$hap_pie_cols, color = NA, show.legend = FALSE) +
  scale_fill_manual(values = hap_assign_col, name = "Haplotype", guide = guide_legend(order = 2)) +
  #scale_fill_manual(values = pop_assign_col, name = "Population") +
  coord_sf(crs = 4326, xlim = c(119, 124.5), ylim = c(20.0, 26.0), expand = FALSE) + # coord_sf(crs = 4326, xlim = c(116.4, 125), ylim = c(20, 26.8), expand = FALSE) +
  
  # === elevation axis custom elements ===
  # outer vertical grid line for the elevation profile
  geom_segment(data = data.frame(x = lwr_bound_axis),
               aes(x = x, y = y_axis_line + y_offset, yend = yend),
               inherit.aes = FALSE, color = "#000000", linewidth = 0.4) +
  # vertical grid lines for the elevation profile
  geom_segment(data = data.frame(x = elev_pos),
               aes(x = x, xend = x, y = y_axis_line + y_offset, yend = yend),
               inherit.aes = FALSE, color = "#9e9e9e", linewidth = 0.4, alpha = 0.3) +
  geom_path(data = elvMax_scaled, aes(x = fake_long, y = latitude, group = 1),
            color = "#3f3f3f", linewidth = 0.4) +
  # axis line spanning only the fake_long strip
  geom_segment(aes(x = lwr_bound_axis, xend = fake_long_max,
                   y = y_axis_line + y_offset, yend = y_axis_line + y_offset),
               linewidth = 0.4) +
  # ticks
  geom_segment(data = data.frame(x = elev_pos),
               aes(x = x, xend = x,
                   y = y_axis_line + y_offset,
                   yend = y_axis_line + y_offset - 0.05),
               inherit.aes = FALSE, linewidth = 0.4) +
  # labels
  geom_text(data = data.frame(x = elev_pos, lab = elev_breaks),
            aes(x = x, y = y_axis_line + y_offset - 0.10, label = lab),
            inherit.aes = FALSE, size = 2.6) +
  # axis title
  annotate("text",
           x = mean(c(lwr_bound_axis, fake_long_max)),
           y = y_axis_line + y_offset - 0.25,
           label = "Elevation (m)",
           size = 3.2, fontface = "bold") +
  
  # === Haplotype network layers ===
  # curve line plot for haplotype links
  geom_curve(data = subset(edges_df, !(x == xend & y == yend)),
             aes(x = x, y = y, xend = xend, yend = yend, linewidth = 0.6),
             curvature = 0.5, colour = "#000000") + #, alpha = 0.6) +
  # straight line plot for haplotype links
  #geom_segment(data = subset(edges_df, !(x == xend & y == yend)),
  #             aes(x = x, y = y, xend = xend, yend = yend, linewidth = 0.6),
  #             colour = "#000000") + #, lineend = "round", alpha = 0.5) +
  scale_linewidth_continuous(range = c(0.2, 1.2), guide = "none") +
  # labels for haplotype step(s)
  geom_label(data = subset(midpoints, !(x == xend & y == yend)),
            aes(x = mx, y = my, label = weight),
            colour = "#000000", fill = "white",label.size = 0.2, size = 3.88) +
  # pie plots for haplotype composition
  geom_scatterpie(data = hap_points_reshaped, # hap_coords
                  aes(x = Longitude, y = Latitude, r = 0.08), 
                  cols = nodes_info$hap_pie_cols, color = "#000000") +
  # labels for haplotypes
  #ggrepel::geom_label_repel(
  #  data = hap_points, #[!(hap_points$Population %in% 
  #                      #hap_points$Population[duplicated(hap_points$Population)]), ],
  #  aes(x = Longitude, y = Latitude, label = Haplotype, fill = Haplotype), # fill = Population
  #  colour = "#000000", label.size = 0.2, show.legend = FALSE,
  #  max.overlaps = Inf,          # allow all haplotypes to be shown
  #  box.padding = 0.25,          # space around label
  #  point.padding = 0.1,         # space around anchor point
  #  min.segment.length = 0,      # always draw leader lines
  #  segment.color = "#000000",
  #  segment.size = 0.4, 
  #  seed = 123                   # reproducible positioning
  #) +
  # === axes, theme, scale bar ===
  scale_x_continuous(limits = c(119, 125), breaks = c(120, 121, 122)) + # scale_x_continuous(limits = c(116.4, 125), breaks = c(118, 120, 122)) +
  scale_y_continuous(limits = c(19.0, 26.0), breaks = c(21, 22, 23, 24, 25)) + # scale_y_continuous(limits = c(20.3, 26.5), breaks = c(21, 22, 23, 24, 25, 26)) +
  annotation_scale(location = "br", width_hint = 0.3, text_cex = 0.8, line_width = 0.4) +
  theme(
    axis.text = element_text(color = "#000000", face = "bold", size = 12),
    axis.title = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "inside",
    legend.position.inside = c(0.00, 1.00),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'transparent'),
    panel.border = element_rect(fill = 'transparent'),
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = 'transparent', color = NA)
  )

#############################################
# Another Plot (when haplotypes overlapped) #
#############################################

hap_order <- c("H4", "H3")
sub_haplo <- ggplot() + 
  # Edges (haplotype links) 
  geom_segment(
    data = subset(edges_df, (x == xend & y == yend)),
    aes(x = factor(from, levels = hap_order), 
        xend = factor(to, levels = hap_order), 
        y = 0, yend = 0
        ), 
    linewidth = 0.6, 
    color = "#000000" 
    ) + 
  scale_linewidth_continuous(range = c(0.2, 1.2), guide = "none") + 
  # Mutation step labels 
  geom_label(
    data = subset(edges_df, (x == xend & y == yend)),
    aes(x = (as.numeric(factor(from, levels = hap_order)) 
             + as.numeric(factor(to, levels = hap_order))) / 2,
        y = 0, label = weight
        ), 
    fill = "#FFFFFF", 
    label.size = 0.2, 
    size = 3.88 
    ) + 
  geom_point(
    data = hap_points[hap_points$Population %in% 
                        hap_points$Population[duplicated(hap_points$Population)], ],
    aes(x = Haplotype, y = 0, fill = Haplotype),
    shape = 21,              
    size = 6.4,              
    color = "#000000",       
    show.legend = FALSE
  ) +
  scale_fill_manual(values = hap_assign_col, name = "Haplotype") +
  scale_x_discrete(drop = FALSE, expand = expansion(mult = c(2.00, 2.00))) + #c(2.00, 2.00))
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    axis.title = element_blank(), 
    panel.background = element_blank(), 
    panel.border = element_blank(),
    panel.grid = element_blank(), 
    plot.background = element_blank() 
    )

####################
# Plot Combination #
####################

main_haplo +
  annotation_custom(
    grob = ggplotGrob(sub_haplo),
    xmin = 119.1, xmax = 120.1,   # position of inset box inside the map (tweak)
    ymin = 20.1, ymax = 21.1      # adjust for where you want it placed
  ) +
  annotate("rect", xmin = 119.1, xmax = 120.1, ymin = 20.1, ymax = 21.1,
           colour = "black", fill = NA, linewidth = 0.4
  )
