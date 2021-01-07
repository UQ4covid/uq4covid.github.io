## load libraries
library(tidyverse)
library(lubridate)
library(sf)
library(areal)
library(patchwork)
library(viridis)

## load shapefiles
ward11 <- st_read("Wards__December_2011__Boundaries_EW_BFC-shp/Wards__December_2011__Boundaries_EW_BFC.shp")
ward19 <- st_read("Wards__December_2019__Boundaries_EW_BFC-shp/Wards__December_2019__Boundaries_EW_BFC.shp")

## make shapefiles valid for matching (solves error if not done)
ward19 <- st_make_valid(ward19)
ward11 <- st_make_valid(ward11)

## set paths
pathToMetaWardsData <- "../../../../../MetaWardsData"
path <- paste0(pathToMetaWardsData, "/model_data/2011Data/")
pathToNewData <- "2011to2019Data"

## load WorkSize tables
WorkSize <- read_table2(paste0(path, "/WorkSize.dat"), col_names = FALSE)
WorkSize19 <- read_table2(paste0(pathToNewData, "/WorkSize19.dat"), col_names = FALSE)

## load lookup tables
Ward_Lookup <- read_csv(paste0(path, "/Ward_Lookup.csv"), col_names = TRUE, guess_max = 3000)
Ward19_Lookup <- read_csv(paste0(pathToNewData, "/Ward19_Lookup.csv"), col_names = TRUE, guess_max = 3000)

## bind shapefiles to lookups
ward19 <- inner_join(ward19, Ward19_Lookup, by = c("wd19cd" = "WD19CD", "wd19nm" = "WD19NM"))
ward11 <- inner_join(ward11, Ward_Lookup, by = c("wd11cd" = "WD11CD", "wd11nm" = "WD11NM"))
ward19 <- inner_join(ward19, WorkSize19, by = c("FID" = "X1"))
ward11 <- inner_join(ward11, WorkSize, by = c("FID" = "X1"))

## produce plots
LADs <- c("City of London", "Exeter", "Hartlepool")
p <- list(NULL)
j <- 1
for(i in 1:length(LADs)) {
    p[[j]] <- filter(ward11, LAD11NM == LADs[i]) %>%
        ggplot() + geom_sf(aes(fill = X2)) + 
        scale_fill_viridis(
            #limits = range(c(ward11$X2[ward11$LAD11NM == LADs[i]], ward19$X2[ward19$LAD19NM == LADs[i]]))
        ) + ggtitle(LADs[i])
    j <- j + 1
    p[[j]] <- filter(ward19, LAD19NM == LADs[i]) %>%
        ggplot() + geom_sf(aes(fill = X2)) + 
        scale_fill_viridis(
            #limits = range(c(ward11$X2[ward11$LAD11NM == LADs[i]], ward19$X2[ward19$LAD19NM == LADs[i]]))
        ) + ggtitle(LADs[i])
    j <- j + 1
}
wrap_plots(p, ncol = 2)
