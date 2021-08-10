## load libraries
library(tidyverse)
library(lubridate)
library(sf)
library(patchwork)
library(viridis)

## load shapefiles
ward11 <- st_read("Wards__December_2011__Boundaries_EW_BFC-shp/Wards__December_2011__Boundaries_EW_BFC.shp")
ward19 <- st_read("Wards__December_2019__Boundaries_EW_BFC-shp/Wards__December_2019__Boundaries_EW_BFC.shp")

## make shapefiles valid for matching (solves error if not done)
ward19 <- st_make_valid(ward19)
ward11 <- st_make_valid(ward11)

## set paths
pathToMetaWardsData <- "../../../../../../MetaWardsData"
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

## some general checks at the LAD level

## extract matching LADs
temp19 <- inner_join(Ward19_Lookup, WorkSize19, by = c("FID" = "X1")) %>%
    group_by(LAD19NM, LAD19CD) %>%
    summarise(WK19 = sum(X2), .groups = "drop")
temp11 <- inner_join(Ward_Lookup, WorkSize, by = c("FID" = "X1")) %>%
    group_by(LAD11NM, LAD11CD) %>%
    summarise(WK11 = sum(X2), .groups = "drop")
temp <- inner_join(temp19, temp11, by = c("LAD19NM" = "LAD11NM")) %>%
    mutate(propdiff = WK11 / WK19) %>%
    mutate(propdiff = ifelse(propdiff < 1, 1 / propdiff, propdiff)) %>%
    arrange(-propdiff)

## produce a few comparison plots
LADs <- temp$LAD19NM[1:3]
p <- list(NULL)
j <- 1
for(i in 1:length(LADs)) {
    p[[j]] <- filter(ward11, LAD11NM == LADs[i]) %>%
        ggplot() + geom_sf(aes(fill = X2)) + 
        scale_fill_viridis(
            limits = range(c(ward11$X2[ward11$LAD11NM == LADs[i]], ward19$X2[ward19$LAD19NM == LADs[i]]))
        ) + ggtitle(paste0(LADs[i], " - 2011")) + labs(fill = "WorkSize")
    j <- j + 1
    p[[j]] <- filter(ward19, LAD19NM == LADs[i]) %>%
        ggplot() + geom_sf(aes(fill = X2)) + 
        scale_fill_viridis(
            limits = range(c(ward11$X2[ward11$LAD11NM == LADs[i]], ward19$X2[ward19$LAD19NM == LADs[i]]))
        ) + ggtitle(paste0(LADs[i], " - 2019")) + labs(fill = "WorkSize")
    j <- j + 1
}
p <- wrap_plots(p, ncol = 2)
ggsave("LADcomps.pdf", p)

## read in LAD shapefiles
lad11 <- st_read("Local_Authority_Districts_(December_2011)_Boundaries_EW_BFC/Local_Authority_Districts_(December_2011)_Boundaries_EW_BFC.shp")
lad19 <- st_read("Local_Authority_Districts_(December_2019)_Boundaries_UK_BFC/Local_Authority_Districts_(December_2019)_Boundaries_UK_BFC.shp")

## check shapefiles match
anti_join(st_drop_geometry(ward11), st_drop_geometry(lad11), by = c("LAD11CD" = "lad11cd"))
anti_join(st_drop_geometry(ward19), st_drop_geometry(lad19), by = c("LAD19CD" = "lad19cd"))

## take a look at LAD-level working populations
temp11 <- st_drop_geometry(ward11) %>%
    group_by(LAD11CD) %>%
    summarise(X2 = sum(X2)) %>%
    {inner_join(lad11, ., by = c("lad11cd" = "LAD11CD"))}
temp19 <- st_drop_geometry(ward19) %>%
    group_by(LAD19CD) %>%
    summarise(X2 = sum(X2)) %>%
    {inner_join(lad19, ., by = c("lad19cd" = "LAD19CD"))}
p1 <- ggplot(temp11) + geom_sf(aes(fill = X2), colour = NA) +
    scale_fill_viridis(
        limits = range(c(temp11$X2, temp19$X2))
    ) + ggtitle("2011") + labs(fill = "WorkSize")   
p2 <- ggplot(temp19) + geom_sf(aes(fill = X2), colour = NA) +
    scale_fill_viridis(
        limits = range(c(temp11$X2, temp19$X2))
    ) + ggtitle("2019") + labs(fill = "WorkSize")
p <- p1 + p2 + plot_layout(guides = "collect")
ggsave("LADWorksize.png", p, width = 10, height = 5)

## load movement data
EW11 <- read_table2(paste0(path, "/EW1.dat"), col_names = FALSE)
EW19 <- read_table2(paste0(pathToNewData, "/EW19.dat"), col_names = FALSE)

## aggregate movements to LADs
temp11 <- st_drop_geometry(ward11) %>%
    select(FID, LAD11CD, LAD11NM)
temp11 <- left_join(EW11, temp11, by = c("X1" = "FID")) %>%
    left_join(temp11, by = c("X2" = "FID")) %>%
    group_by(LAD11NM.x, LAD11NM.y) %>%
    summarise(X3 = sum(X3), .groups = "drop") %>%
    mutate(X3 = log(X3))
temp19 <- st_drop_geometry(ward19) %>%
    select(FID, LAD19CD, LAD19NM)
temp19 <- left_join(EW19, temp19, by = c("X1" = "FID")) %>%
    left_join(temp19, by = c("X2" = "FID")) %>%
    group_by(LAD19NM.x, LAD19NM.y) %>%
    summarise(X3 = sum(X3), .groups = "drop") %>%
    mutate(X3 = log(X3))

## produce a few comparison plots
LADs <- temp$LAD19NM[1:3]
p <- list(NULL)
j <- 1
for(i in 1:length(LADs)) {
    p1 <- filter(temp11, LAD11NM.x == LADs[i]) %>%
        select(-LAD11NM.x) %>%
        {inner_join(lad11, ., by = c("lad11nm" = "LAD11NM.y"))} 
    p2 <- filter(temp19, LAD19NM.x == LADs[i]) %>%
        select(-LAD19NM.x) %>%
        {inner_join(lad19, ., by = c("lad19nm" = "LAD19NM.y"))}
    p1 <- ggplot(p1) + geom_sf(aes(fill = X3)) + 
        scale_fill_viridis(
            limits = range(c(p1$X3, p2$X3))
        ) + ggtitle(paste0("From ", LADs[i], " - 2011")) + labs(fill = "Log-movements")
    p2 <- ggplot(p2) + geom_sf(aes(fill = X3)) + 
        scale_fill_viridis(
            limits = range(c(p1$X3, p2$X3))
        ) + ggtitle(paste0("From ", LADs[i], " - 2019")) + labs(fill = "Log-movements")
    p[[i]] <- p1 + p2 + plot_layout(guides = "collect")
}
p <- wrap_plots(p, ncol = 1)
ggsave("LADmovements.png", p, width = 10, height = 15)

