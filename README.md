# distribution_map
plot sampling sites on a map and an elevation profile

Follow the script, prepare your input files, and then keep pressing run :)

# **----Link of the TIFF file----**
https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1
<br.>You can manually enter the coordinates according to your specific needs.
<br.>The bound (116.5, 20.5, 122.2, 26.4) includes Dongsha Island, the Matsu Islands, the Kinmen Islands, Taiwan Island, Green Island, and the Orchid Islands.

# **----Things to keep in mind----**

The plot may vary for different users; you may need to adjust it accordingly.

Especially:
##### outer vertical grid line for the elevation profile
geom_segment(data = data.frame(x = lwr_bound_axis), aes(x = x, y = y_axis_line + y_offset, yend = 26.5), inherit.aes = FALSE, color = "#000000", linewidth = 0.4) + 
##### vertical grid lines for the elevation profile
geom_segment(data = data.frame(x = elev_pos), aes(x = x, xend = x, y = y_axis_line + y_offset, yend = 26.5), inherit.aes = FALSE, color = "#9e9e9e", linewidth = 0.4, alpha = 0.3) + 
##### line plot for the elevation profile
geom_path(data = elvMax_scaled, aes(x = fake_long, y = latitude, group = 1), color = "#3f3f3f", linewidth = 0.4) + 
##### axis line spanning only the fake_long strip
geom_segment(aes(x = lwr_bound_axis, xend = fake_long_max, y = y_axis_line + y_offset, yend = y_axis_line + y_offset), linewidth = 0.4) +
##### axis ticks
geom_segment(data = data.frame(x = elev_pos), aes(x = x, xend = x, y = y_axis_line + y_offset, yend = y_axis_line + y_offset - 0.05), inherit.aes = FALSE, linewidth = 0.4) +
##### labels
geom_text(data = data.frame(x = elev_pos, lab = elev_breaks), aes(x = x, y = y_axis_line + y_offset - 0.10, label = lab), inherit.aes = FALSE, size = 2.5) + 
##### axis title
annotate("text", x = mean(c(lwr_bound_axis, fake_long_max)), y = y_axis_line + y_offset - 0.25, label = "Elevation (m)", size = 3.2, fontface = "bold")

This section is for creating the secondary x-axis for the elevation profile.

# **----========Section Break========----**

For the theme, you can modify it based on your preferences and also your computer resolution.
  theme(
    axis.text = element_text(color = "#000000", face = "bold"),
    axis.title = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "inside",
    legend.position.inside = c(0.00, 1.00),
    legend.justification = c("left", "top"),
    legend.box.just = "right",
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'transparent'),
    panel.border = element_rect(fill = 'transparent'),
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = 'transparent', color = NA)
  )

  
It is easy to follow, I think xD
