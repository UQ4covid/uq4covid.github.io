## load libraries
library(tidyverse)
library(sf)

## load shapefile
lad19 <- st_read("Local_Authority_Districts_(December_2019)_Boundaries_UK_BUC/Local_Authority_Districts_(December_2019)_Boundaries_UK_BUC.shp")

## check shapefile has single entry per LAD CD and NM
stopifnot(
    group_by(lad19, lad19cd, lad19nm) %>%
        count() %>%
        pluck("n") %>%
        {all(. == 1)}
)

## check FIDs complete
stopifnot(all((lad19$objectid - 1:nrow(lad19)) == 0))

## set paths to ward lookup
path <- paste0("../../../../../../MetaWardsData/model_data/2011to2019Data/")

## load ward lookup in data repo
wards <- read_csv(paste0(path, "Ward19_Lookup.csv"))

## check for any missing wards in new lookup
anti_join(wards, lad19, by = c("LAD19CD" = "lad19cd", "LAD19NM" = "lad19nm"))

## remove extraneous LADs from shapefile
lad19 <- semi_join(lad19, wards, by = c("lad19cd" = "LAD19CD", "lad19nm" = "LAD19NM")) %>%
    mutate(objectid = 1:n())

## add LAD FID
wards <- left_join(wards, select(st_drop_geometry(lad19), objectid, lad19cd, lad19nm), 
    by = c("LAD19CD" = "lad19cd", "LAD19NM" = "lad19nm"))

## load different data sets
CBB2019 <- read_table2(paste0(path, "CBB2019.dat"), col_names = FALSE)
EW19 <- read_table2(paste0(path, "EW19.dat"), col_names = FALSE)
PlayMatrix19 <- read_table2(paste0(path, "PlayMatrix19.dat"), col_names = FALSE)
PlaySize19 <- read_table2(paste0(path, "PlaySize19.dat"), col_names = FALSE)
WorkSize19 <- read_table2(paste0(path, "WorkSize19.dat"), col_names = FALSE)

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## PlayMatrix19 is the proportion of players (X3) moving from ward X1 to ward X2 in play movements

## check proportions sum to one
group_by(PlayMatrix19, X1) %>%
    summarise(X3 = sum(X3), .groups = "drop") %>%
    pluck("X3") %>%
    summary()

## calculate average number of movements whilst retaining population size
PlayMatrix19 <- left_join(PlayMatrix19, rename(PlaySize19, n = X2), by = "X1") %>%
    mutate(np = n * X3) %>%
    group_by(X1) %>%
    mutate(np = smart_round(np))

## check sums match
select(PlayMatrix19, X1, n, np) %>%
    group_by(X1) %>%
    summarise(n = mean(n), np = sum(np), .groups = "drop") %>%
    mutate(diff = n - np) %>%
    pluck("diff") %>%
    summary()

## now merge movements to LAD level and recalculate proportions
PlayMatrix19 <- select(PlayMatrix19, X1, X2, np) %>%
    left_join(select(wards, FID, objectid), by = c("X1" = "FID")) %>%
    mutate(X1 = objectid) %>%
    select(-objectid) %>%
    left_join(select(wards, FID, objectid), by = c("X2" = "FID")) %>%
    mutate(X2 = objectid) %>%
    select(-objectid) %>%
    group_by(X1, X2) %>%
    summarise(np = sum(np), .groups = "drop") %>%
    group_by(X1) %>%
    mutate(X3 = np / sum(np)) %>%
    ungroup() %>%
    select(!np) %>%
    arrange(X1, X2)

## check proportions match
group_by(PlayMatrix19, X1) %>%
    summarise(X3 = sum(X3), .groups = "drop") %>%
    pluck("X3") %>%
    summary()

## EW19 is the number of individuals (X3) who move from ward A (X1) to ward B (X2)
EW19 <- left_join(EW19, select(wards, FID, objectid), by = c("X1" = "FID")) %>%
    mutate(X1 = objectid) %>%
    select(-objectid) %>%
    left_join(select(wards, FID, objectid), by = c("X2" = "FID")) %>%
    mutate(X2 = objectid) %>%
    select(-objectid) %>%
    group_by(X1, X2) %>%
    summarise(X3 = sum(X3), .groups = "drop") %>%
    arrange(X1, X2)

## WorkSize19 is the number of workers (X2) in ward (X1)
WorkSize19 <- left_join(WorkSize19, select(wards, FID, objectid), by = c("X1" = "FID")) %>%
    mutate(X1 = objectid) %>%
    select(-objectid) %>%
    group_by(X1) %>%
    summarise(X2 = sum(X2), .groups = "drop") %>%
    arrange(X1)

## PlaySize19 is the number of players (X2) in ward (X1)
PlaySize19 <- left_join(PlaySize19, select(wards, FID, objectid), by = c("X1" = "FID")) %>%
    mutate(X1 = objectid) %>%
    select(-objectid) %>%
    group_by(X1) %>%
    summarise(X2 = sum(X2), .groups = "drop") %>%
    arrange(X1)

## CBB2019 is point locations of wards, so now convert to LADs
CBB2019 <- st_centroid(lad19) %>% 
    {do.call("rbind", st_geometry(.))} %>%
    {cbind(lad19$objectid, .)} %>%
    as_tibble() %>%
    arrange(V1)
colnames(CBB2019) <- paste0("X", 1:ncol(CBB2019))

## create directory to write out data
pathToNewData <- "2019LADData"
dir.create(pathToNewData)

## create description file
description <- c(
    "{ \"name\"               : \"2019LADData\",",
    " \"version\"            : \"Jan 25 2022\",",
    " \"author(s)\"          : \"TJ McKinley\",",
    "  \"contact(s)\"         : \"t.mckinley@exeter.ac.uk\",",
    "  \"work\"               : \"EW19.dat\",",
    "  \"work_size\"          : \"WorkSize19.dat\",",
    "  \"play\"               : \"PlayMatrix19.dat\",",
    "  \"play_size\"          : \"PlaySize19.dat\",",
    "  \"position\"           : \"CBB2019.dat\",",
    "  \"coordinates\"        : \"x/y\",",
    "  \"coordinate_units\"   : \"m\",",
    "  \"lookup\"             : \"LAD19_Lookup.csv\",",
    "  \"lookup_columns\"     : {\"code\":1, \"name\":2}",
    "}")
writeLines(description, paste0(pathToNewData, "/description.json"))

## write all other files
write_delim(EW19, paste0(pathToNewData, "/EW19.dat"), col_names = FALSE)
write_delim(WorkSize19, paste0(pathToNewData, "/WorkSize19.dat"), col_names = FALSE)
write_delim(PlayMatrix19, paste0(pathToNewData, "/PlayMatrix19.dat"), col_names = FALSE)
write_delim(PlaySize19, paste0(pathToNewData, "/PlaySize19.dat"), col_names = FALSE)
write_delim(CBB2019, paste0(pathToNewData, "/CBB2019.dat"), col_names = FALSE)
write_csv(select(st_drop_geometry(lad19), FID = objectid, LAD19CD = lad19cd, LAD19NM = lad19nm), paste0(pathToNewData, "/LAD19_Lookup.csv"))
write_csv(select(st_drop_geometry(lad19), FID = objectid, LAD19CD = lad19cd, LAD19NM = lad19nm), "../../inputs/LAD19_Lookup.csv")

## write new shapefile
dir.create("LAD19_shapefile")
st_write(lad19, "LAD19_shapefile/LAD19_shapefile.shp")
