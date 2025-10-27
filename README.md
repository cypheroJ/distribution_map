# hap_map.R
Plot haplotypes on a map and an elevation profile

<img width="2940" height="1840" alt="image" src="https://github.com/user-attachments/assets/a55527bf-cacb-4fde-af6e-96765ce5d208" />


Files provided for instructions:
1. haplo_co1.fasta
2. haplo_pop.csv
3. samp_site.R
You can follow the files that I provided and modify them for your own purposes.


# **----Link of the shp files----**
https://data.moi.gov.tw/MoiOD/System/DownloadFile.aspx?DATA=72874C55-884D-4CEA-B7D6-F60B0BE85AB0 
This link provides the country border of Taiwan, including islands from South China Sea, Kinmen, Matsu (Lienchiang County), Penghu, Orchid Islands, Green Island, and the main island.


# **----Link of the TIFF file----**
https://portal.opentopography.org/raster?opentopoID=OTSRTM.082015.4326.1  
You can manually enter the coordinates according to your specific needs.  
The bound (116.5, 20.5, 122.2, 26.4) includes Dongsha Island, the Matsu Islands, the Kinmen Islands, Taiwan Island, Green Island, and the Orchid Islands.
We will need two bounds if we have the data from the South China Sea islands.
(116.5, 21.45014, 124.561, 26.38542) & (114.359, 10.37153, 116.9999, 21.45014)


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

# **Resolution selection for map** #
geom_spatraster(data = elv_masked, maxcell = 3.8e7)  
You can try geom_spatraster(data = elv_masked) <SpatRaster> resampled to 501120 cells.  
maxcell = 3.8e7 <SpatRaster> resampled to 38008320 cells.  
or maxcell = 1e6 or lower or higher (not more than 3.8e7).  

# **========Section Break========**

For the scalebar, you can modify it by changing the location:
"br" bottom-right, "bl" bottom-left, "tr" top-right, "tl" top-left
annotation_scale(location = "br", width_hint = 0.3, text_cex = 0.8, line_width = 0.4)

For the theme, you can modify it based on your preferences and computer resolution.  
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
