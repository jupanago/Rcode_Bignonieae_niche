## This script is for producing the environment of AoEs from AoEs' polygons

# NOTE: This should not be run again. Data from environment extraction in:
# "./data/biogeography/AoE_raster/AoE_Chelsa10m_bio.csv"

# ### Rasterizing
# 
# ### Function ndm_rasterize
# ##
# ## Arguments: 
# ## (i) 'filepath': filepath of the NDM grid output.
# ## (ii) 'dg': spatial scale in degrees. This argument goes in the argument 'width' of the function rgeos::gBuffer, which is then divided by 2. 
# ## (iii) 'save_area.shp': default = FALSE. Save area as shapefile using function rgdal::writeOGR.
# ## (iv) 'save.directory':  
# ## (v) area.shp.name: default = "Area_id". Name of the area shapefile. 
# ## (vi) env_ras: path to raster file that will be used as template
# ## (vii) region_shp_dns: path to shapefile of the geographical area in which the AoE occurs
# ## (viii) layer_name: name of the geographical area in which the AoE occurs
# 
# ndm_rasterize <- function (filepath, degree, save.area.shp = FALSE, save.raster = FALSE, 
#                            save.directory, area.shp.name = "Area_id", env_ras, region_shp_dns, layer_name) {
#   
#   require(raster)
#   require(sp)
#   require(rgdal)
#   require(rgeos)
#   require(sf)
#   require(rmapshaper)
#   
#   ### Loading data
#   
#   ## Loading environmental template raster
#   ras_env <- raster(env_ras) 
#   ## Loading shapefile to define the geographical range
#   continent <- readOGR(dsn = region_shp_dns, layer = layer_name)
#   
#   ########### Parsing grid.txt file from NDM ########### 
#   
#   ## x is a vector whose elements correspond to each line in the gridxy.txt archive.
#   x <- readLines(con = filepath)
#   
#   ## Initializing data.frame
#   ## ini_df colums: "id" = index for coordinates belongind to the same grid, "ord" = points joining order, 
#   ## "long" = longitude, "lat" = latitude. 
#   ## NOTE: "ord" can be deleted. It was created to account for joining order in drawing the cell polygons if necessary.
#   
#   ini_df <- data.frame(id = character(),
#                        ord = double(),
#                        long = character(),
#                        lat = character())
#   
#   ## Counting for n = id of the points belonging to the same cell, and y = order for joining points in a cell.
#   n <- 0
#   y <- 1
#   
#   ## Loop for organizing vector elements in their corresponding places inside the data.frame 
#   ## Regex: "(^[0-9]+) (.[0-9.]+) (.[0-9.]+)" describes the pattern of element of the vector that
#   ## contain 3 numbers, one identifier and two coordinates.
#   ## Regex: "^(.[0-9.]+) (.[0-9.]+)" describes the patter of the vector element that contains only
#   ## both coordinates.
#   
#   for (i in seq_along(x)) {
#     
#     if (grepl("(^[0-9]+) (.[0-9.]+) (.[0-9.]+)", x[i])) {
#       ini_df <- rbind(ini_df, data.frame(id = sub("(^[0-9]+) (.[0-9.]+) (.[0-9.]+)", "\\1", x[i]),
#                                          ord = y,
#                                          long = sub("(^[0-9]+) (.[0-9.]+) (.[0-9.]+)", "\\2", x[i]),
#                                          lat = sub("(^[0-9]+) (.[0-9.]+) (.[0-9.]+)", "\\3", x[i])))
#       
#       n <- n + 1
#       y <- y + 1
#       
#     }
#     
#     if (grepl("^(.[0-9.]+) (.[0-9.]+)$", x[i])) {
#       ini_df <- rbind(ini_df, data.frame(id = n,
#                                          ord = y,
#                                          long = sub("^(.[0-9.]+) (.[0-9.]+)$", "\\1", x[i]),
#                                          lat = sub("^(.[0-9.]+) (.[0-9.]+)$", "\\2", x[i])))
#       y <- y + 1
#       
#     }
#     
#     # Condition to restart the count of y for each cell.
#     if (y == 6) {
#       
#       y
#       y <- 1
#       
#     }
#     
#   }
#   
#   ## Changing columns class from character to numeric. gsub, grep, and grepl work with characters. 
#   ini_df[,1] <- as.double(as.character(ini_df[,1]))
#   ini_df[,3] <- as.double(as.character(ini_df[,3]))
#   ini_df[,4] <- as.double(as.character(ini_df[,4]))
#   
#   ini_df
#   
#   ########### Creating a grid shapefile ########### 
#   
#   centroid_df <- data.frame(id = double(),
#                             long = double(),
#                             lat = double())
#   
#   for (i in unique(ini_df$id)) {
#     
#     int_df <- subset(ini_df, subset = id == i)
#     centroid_df <- rbind(centroid_df, data.frame(id = i,
#                                                  long = mean(int_df[1:4, 3]),
#                                                  lat = mean(int_df[1:4, 4])))
#     
#     centroid_df
#     
#   }
#   
#   df = data.frame(id = centroid_df[,1])
#   centroid_spdf <- sp::SpatialPointsDataFrame(centroid_df[,2:3], df, proj4string = CRS(as.character(NA)))
#   
#   ## Creating the grid
#   grid <- rgeos::gBuffer(centroid_spdf, byid = T, id = centroid_df[,1],
#                          width = degree/2, capStyle = "SQUARE")
#   
#   
#   if (save.area.shp == TRUE) { # Saving area shapefile.
#     
#     rgdal::writeOGR(grid, dsn = save.directory, layer = area.shp.name, driver = "ESRI Shapefile")
#     
#   }
#   
#   
#   ########## Clipping the geographical area to the extent of the environmental raster using a bounding box ##########
#   
#   ## Function gClip() from https://www.r-bloggers.com/clipping-spatial-data-in-r/ 
#   gClip <- function(shp, bb){  
#     if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
#     else b_poly <- as(extent(bb), "SpatialPolygons")
#     gIntersection(shp, b_poly, byid = T)
#   }
#   ## Clipping the shapefile
#   continent_clip <- gClip(shp = continent, bb = bbox(ras_env))
#   
#   ### Cliping grid to the geographical area
#   continent_clip <- st_as_sf(continent_clip) 
#   st_crs(continent_clip) <- "+proj=longlat +datum=WGS84 +no_defs"
#   
#   grid_ds <- grid %>% ms_dissolve() %>% st_as_sf() 
#   st_crs(grid_ds) <- "+proj=longlat +datum=WGS84 +no_defs"
#   
#   grid <- st_intersection(grid_ds, continent_clip )
#   
#   ######### Rasterizing AoEs #########
#   
#   ### Rasterizing: Creating raster template for the sea
#   
#   ## Using ras_env as a template
#   raster_sea <- as(raster::extent(ras_env), "SpatialPolygons") 
#   ## Defining projection
#   proj4string(raster_sea) <- crs(st_crs(ras_env)$proj4string) 
#   ## Defining resolution
#   raster_sea  <- raster(raster_sea, resolution = res(ras_env))
#   ## Setting values to NA. This will be the value for the sea
#   raster_sea[] <- NA 
#   
#   ### Rasterizing: continent_clip using ras_env as template
#   continent_ras <- rasterize(continent_clip, ras_env, field = 0) # Continent values at 0
#   
#   ### Rasterizing: grid of area of endemism
#   grid_ras <- rasterize(grid, continent_ras, field = 1)
#   
#   ### Merging rasters into one layer
#   AoE_raster <- merge(grid_ras, continent_ras, raster_sea, overwrite = FALSE)
#   
#   # Confirming everything is OK!
#   # plot(AoE_raster)
#   # unique(values(AoE))
#   # res(AoE) == res(ras_env)
#   # origin(AoE) == origin(ras_env)
#   # st_crs(AoE) == st_crs(ras_env)
#   
#   ## Saving raster
#   
#   if (save.raster == TRUE) {
#     
#     raster::writeRaster(AoE_raster, filename = paste(save.directory, area.shp.name, sep = "") , format = "GTiff", overwrite = TRUE)
#     
#   }
#   
#   
#   
#   
# }
# 
# 
# # ndm_rasterize(filepath = "./data/GIS/40perc/2dgD_CA17_NWCol-grid.txt", degree = 2, save.raster = TRUE, 
# #               save.directory = "./output/AoE_raster/", area.shp.name = "perfect", 
# #               env_ras = "./data/env/chelsa10m/bio/bio_1.tif", region_shp_dns = "./data/GIS/Amer_land.shp",
# #               layer_name = "Amer_land")
# 
# # All areas at once.
# 
# ndm_listgrid.files <- function (path, degree) {
#   
#   # This function is specifically designed to my own code of file names: 
#   # degreeAnalysis_ConsensusArea#_AreaName-grid.txt => Example: "1dgT_CA0_AF-grid.txt". General Regex = "^.*grid.txt 
#   # Arguments: 
#   # (i) 'path' is the filepath whre the ndm grid outputs are stored. 
#   # (ii) 'degree': It is the spatial degree of analysis which is refered at the beggining of the archive name.
#   
#   if (degree == 1) {
#     
#     files1dg <<- list.files(path = path, pattern = "^1.*grid.txt", full.names = TRUE)
#     anames1dg <<- gsub(pattern = "^.*(1dg.*).txt", replacement = "\\1", files1dg)
#     message(paste(length(files1dg), "files listed in 'files1dg'.", length(anames1dg), "area names generated and saved in 'anames1dg'."))
#     
#   } 
#   
#   if (degree == 2) {
#     
#     files2dg <<- list.files(path = path, pattern = "^2.*grid.txt", full.names = TRUE)
#     anames2dg <<- gsub(pattern = "^.*(2dg.*).txt", replacement = "\\1", files2dg)
#     message(paste(length(files2dg), "files listed in 'files2dg'.", length(anames2dg), "area names generated and saved in 'anames2dg'."))
#     
#   } 
#   
#   if (degree == 3) {
#     
#     files3dg <<- list.files(path = path, pattern = "^3.*grid.txt", full.names = TRUE)
#     anames3dg <<- gsub(pattern = "^.*(3dg.*).txt", replacement = "\\1", files3dg)
#     message(paste(length(files3dg), "files listed in 'files3dg'.", length(anames3dg), "area names generated and saved in 'anames3dg'."))
#     
#   } 
#   
# }
# ndm_listgrid.files("./data/GIS/40perc/", degree = 2)
# 
# library(purrr)
# 
# #Chelsa10m
# pmap(list(file = files2dg,
#           degree = 2,
#           save.area.shp = FALSE,
#           save.raster = TRUE,
#           save.directory = "./data/AoE_raster/Chelsa10m/",
#           area.shp.name = anames2dg,
#           env_ras = "./data/env/chelsa10m/bio/bio_1.tif",
#           region_shp_dns = "./data/GIS/Amer_land.shp",
#           layer_name = "Amer_land"), ndm_rasterize)
# 
# 
# ##### Extracting environmental datra
# 
# #### Function ndm_extract
# ##
# ## Arguments:
# ##  x: file path to the directory where the AoE rasters are stored
# ##  y: file path to the directory where the environmental predictors are stored
# ##  save_dns: directory to save the dataframe as .cvs
# ##  cvs_name: name of the cvs file
# 
# ndm_envextract <- function(x, y, save_dns, cvs_name = "AoE_env") {
#   
#   require(raster)
#   require(dplyr)
#   
#   ## Load files
#   #AoE
#   files_a <- list.files(x, pattern = "*grid.tif", full.names = TRUE)
#   AoE_raster <- stack(files_a)
#   names(AoE_raster) <- sub("X2dgD_(CA[0-9]+).*", paste("D2", "\\1", sep = ""), names(AoE_raster))
#   
#   #Environmental predictors
#   files_p <- list.files(y, full.names = T, pattern = "*.tif$")
#   predictors <- stack(files_p)
#   
#   #Creating dataframe 
#   col_names <- c("AoE_name", "ID", "cells", names(predictors))
#   ini_df <- data.frame(matrix(ncol = length(col_names), nrow = 0))
#   names(ini_df) <- col_names
#   
#   # Extracting values
#   for (i in 1:length(files_a)) {
#     
#     AoE_points <- xyFromCell(AoE_raster[[i]], which(AoE_raster[[i]][] == 1), spatial = TRUE) #Intermediate step: raster as spatial points
#     AoE_env <- raster::extract(predictors, AoE_points, cellnumber = T, df = T)
#     AoE_env$AoE_name <- names(AoE_raster[[i]])
#     
#     ini_df <- bind_rows(ini_df, AoE_env)
#     
#   }
#   
#   write.csv(ini_df, paste(save_dns, cvs_name, ".csv", sep = ""))
#   
#   return(ini_df)
#   
# }
# 
# # aoe_chelsa10m_bio <- ndm_envextract(x = "./output/AoE_raster/Chelsa10m/",
# #                                     y = "./data/env/chelsa10m/bio/",
# #                                     save_dns = "./output/AoE_raster/Chelsa10m/",
# #                                     cvs_name = "AoE_Chelsa10m_bio")
# 

