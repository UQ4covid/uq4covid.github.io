## load libraries
library(tidyverse)
library(lubridate)
library(sf)
library(areal)

## set paths
pathToMetaWardsData <- "../../../../../../MetaWardsData"
path <- paste0(pathToMetaWardsData, "/model_data/2011Data/")

#############################################################
######  FIND SHAPEFILE INTERSECTIONS AND SHARED AREAS  ######
#############################################################

## load shapefiles
ward11 <- st_read("Wards__December_2011__Boundaries_EW_BFC-shp/Wards__December_2011__Boundaries_EW_BFC.shp")
ward19 <- st_read("Wards__December_2019__Boundaries_EW_BFC-shp/Wards__December_2019__Boundaries_EW_BFC.shp")

## load lookups
Ward19_Lookup <- read_csv("Ward_to_Local_Authority_District_(December_2019)_Lookup_in_the_United_Kingdom.csv")
Ward_Lookup <- read_csv(paste0(path, "Ward_Lookup.csv"), guess_max = 3000) %>%
  arrange(FID)

## check ward lookup has single entry per ward CD and NM
stopifnot(
  group_by(Ward_Lookup, WD11CD, WD11NM) %>%
    count() %>%
    pluck("n") %>%
    {all(. == 1)}
)
## check FIDs complete
stopifnot(all((Ward_Lookup$FID - 1:nrow(Ward_Lookup)) == 0))

## load Ward19 lookup has single entry per ward CD and NM
stopifnot(
  group_by(Ward19_Lookup, WD19CD, WD19NM) %>%
    count() %>%
    pluck("n") %>%
    {all(. == 1)}
)
## check FIDs complete
stopifnot(all((Ward19_Lookup$FID - 1:nrow(Ward19_Lookup)) == 0))

## check names match between shapefiles and lookups
temp <- full_join(ward11, Ward_Lookup, by = c("wd11cd" = "WD11CD", "wd11nm" = "WD11NM"), keep = TRUE)
filter(temp, is.na(wd11cd))
filter(temp, is.na(WD11CD)) 

## these are typos and can be fixed manually
Ward_Lookup <- mutate(Ward_Lookup, WD11NM = ifelse(WD11NM == "Ockbrook and Borrowash", "Ockbrook And Borrowash", WD11NM)) %>%
  mutate(WD11NM = ifelse(WD11NM == "Felin-fâch", "Felin-fƒch", WD11NM))

## check names match between shapefiles and lookups
temp <- full_join(ward11, Ward_Lookup, by = c("wd11cd" = "WD11CD", "wd11nm" = "WD11NM"), keep = TRUE)
stopifnot(nrow(filter(temp, is.na(wd11cd))) == 0)
stopifnot(nrow(filter(temp, is.na(WD11CD))) == 0)

## remove NI and Scottish wards
Ward19_Lookup <- filter(Ward19_Lookup, !str_detect(WD19CD, "^N")) %>%
  filter(!str_detect(WD19CD, "^S")) %>%
  arrange(FID) %>%
  mutate(FID = 1:n())

## check names match between shapefiles and lookups
temp <- full_join(ward19, Ward19_Lookup, by = c("wd19cd" = "WD19CD", "wd19nm" = "WD19NM"), keep = TRUE)
filter(temp, is.na(wd19cd))
filter(temp, is.na(WD19CD))

## these are typos and can be fixed manually
Ward19_Lookup <- mutate(Ward19_Lookup, WD19NM = ifelse(WD19NM == "Canolbarth Môn", "Canolbarth Mwn", WD19NM)) %>%
  mutate(WD19NM = ifelse(WD19NM == "Felin-fâch", "Felin-fach", WD19NM)) %>%
  mutate(WD19NM = ifelse(WD19NM == "Llifôn", "Llifon", WD19NM))

## check names match between shapefiles and lookups
temp <- full_join(ward19, Ward19_Lookup, by = c("wd19cd" = "WD19CD", "wd19nm" = "WD19NM"), keep = TRUE)
stopifnot(nrow(filter(temp, is.na(wd19cd))) == 0)
stopifnot(nrow(filter(temp, is.na(WD19CD))) == 0)

## make shapefiles valid for matching (solves error if not done)
ward19 <- st_make_valid(ward19)
ward11 <- st_make_valid(ward11)

## extract intersection areas and weights, where weights are proportional
## to area of overlap between wards
wardweights <- ward19 %>%
  aw_intersect(source = ward11, areaVar = "area") %>%
  aw_total(source = ward11, id = wd11cd, areaVar = "area", totalVar = "totalArea",
           type = "extensive", weight = "sum") %>%
  aw_weight(areaVar = "area", totalVar = "totalArea", areaWeight = "areaWeight") %>%
  st_drop_geometry()

##################################################
######        READ IN MOVEMENT DATA         ######
##################################################

## load different data sets
CBB2011 <- read_table2(paste0(path, "CBB2011.dat"), col_names = FALSE)
EW1 <- read_table2(paste0(path, "EW1.dat"), col_names = FALSE)
PlayMatrix <- read_table2(paste0(path, "PlayMatrix.dat"), col_names = FALSE)
PlaySize <- read_table2(paste0(path, "PlaySize.dat"), col_names = FALSE)
# seeds <- read_table2(paste0(path, "seeds.dat"), col_names = FALSE) ## not currently used
WorkSize <- read_table2(paste0(path, "WorkSize.dat"), col_names = FALSE)

## CBB2011, Ward_Lookup, WorkSize and PlaySize correspond to a
## ward, indexed by the first column, and some movements of 
## positions in the other columns (I think)

## full details here: https://metawards.org/fileformats/network.html

## check EW1 numbers match to WorkSize
temp <- group_by(EW1, X1) %>%
  summarise(X3 = sum(X3)) %>%
  full_join(WorkSize, by = "X1")
stopifnot(all(temp$X3 == temp$X2))
rm(temp)

## check PlayMatrix proportions add up to one across wards
temp <- PlayMatrix %>%
  group_by(X1) %>%
  summarise(X3 = sum(X3)) %>%
  mutate(eq = map_lgl(X3, ~all.equal(., 1)))
stopifnot(all(temp$eq))
rm(temp)

##################################################
######  AGGREGATE 2011 DATA TO 2019 WARDS   ######
##################################################

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}

## IMPLEMENTATION NOTE: due to rounding errors it is challenging to get all the 
## movements converted whilst still matching across the different sources
## of data - the approach we use here is to work with the movements (EW1 / PlayMatrix)
## and convert those, and then re-construct the WorkSize and PlaySize from those
## to make sure everything matches up in MetaWards. Otherwise, converting e.g. WorkSize
## and then EW1 separately, means that the total counts won't quite match up

## aggregate EW1 (split into multiple pipes because of large memory constraints)

## firstly, map X1 from WD11 to WD19
## (smart_round() ensures total no. of movements match as
## long as grouping is done, even though grouping
## takes a long time)
EW1_new <- EW1 %>%
  left_join(
    select(Ward_Lookup, FID, WD11CD, WD11NM),
    by = c("X1" = "FID")) %>%
  inner_join(
    select(wardweights, wd11cd, wd11nm, wd19cd, wd19nm, areaWeight),
    by = c("WD11CD" = "wd11cd", "WD11NM" = "wd11nm")) %>%
  group_by(X1, X2) %>%
  mutate(X3 = smart_round(X3 * areaWeight)) %>%
  ungroup() %>%
  filter(X3 > 0) %>%
  select(wd19cd, wd19nm, X1, X2, X3)

## check total number of movements out of old X1s match
## after redistributing
stopifnot(group_by(EW1_new, X1) %>% 
  summarise(X3 = sum(X3)) %>%
  inner_join(WorkSize, by = "X1") %>%
  mutate(diff = abs(X3 - X2)) %>%
  pluck("diff") %>%
  {all(. == 0)})

## convert to new X1 based on WD19
EW1_new <- EW1_new %>%
  inner_join(Ward19_Lookup, by = c("wd19cd" = "WD19CD", "wd19nm" = "WD19NM")) %>%
  select(X1 = FID, X2, X3) %>%
  group_by(X1, X2) %>%
  summarise(X3 = sum(X3)) %>%
  ungroup()

## secondly, map X2 from WD11 to WD19
## (smart_round() ensures total no. of movements match as
## long as grouping is done, even though grouping
## takes a long time)
EW1_new <- EW1_new %>%
  left_join(
    select(Ward_Lookup, FID, WD11CD, WD11NM),
    by = c("X2" = "FID")) %>%
  inner_join(
    select(wardweights, wd11cd, wd11nm, wd19cd, wd19nm, areaWeight),
    by = c("WD11CD" = "wd11cd", "WD11NM" = "wd11nm")) %>%
  group_by(X1, X2) %>%
  mutate(X3 = smart_round(X3 * areaWeight)) %>%
  ungroup() %>%
  filter(X3 > 0) %>%
  select(wd19cd, wd19nm, X1, X2, X3)

## check total number of movements into old X2s match
## after redistributing
stopifnot(group_by(EW1_new, X2) %>% 
  summarise(X3 = sum(X3)) %>%
  inner_join(
      EW1 %>%
        group_by(X2) %>%
        summarise(X3 = sum(X3)), 
    by = "X2") %>%
  mutate(diff = abs(X3.x - X3.y)) %>%
  pluck("diff") %>%
  {all(. == 0)})

## convert to new X2 based on WD19
EW1_new <- EW1_new %>%
  inner_join(Ward19_Lookup, by = c("wd19cd" = "WD19CD", "wd19nm" = "WD19NM")) %>%
  select(X1, X2 = FID, X3) %>%
  group_by(X1, X2) %>%
  summarise(X3 = sum(X3)) %>%
  ungroup() %>%
  arrange(X1, X2)

## check total movements match
stopifnot(sum(EW1_new$X3) == sum(EW1$X3))

## now generate new WorkSize
WorkSize_new <- EW1_new %>%
  group_by(X1) %>%
  summarise(X3 = sum(X3)) %>%
  select(X1, X2 = X3) %>%
  arrange(X1)

## create new PlayMatrix (split into multiple pipes because of large memory constraints)

## firstly, create equivalent of EW1 for Players
PW1 <- inner_join(PlayMatrix, PlaySize, by = "X1") %>%
  group_by(X1) %>%
  mutate(X3 = smart_round(X2.y * X3)) %>%
  ungroup() %>%
  select(X1, X2 = X2.x, X3) %>%
  arrange(X1, X2)

## check PW1 numbers match to PlaySize in case of rounding errors
temp <- group_by(PW1, X1) %>%
  summarise(X3 = sum(X3)) %>%
  full_join(PlaySize, by = "X1")
stopifnot(all(temp$X3 == temp$X2))
rm(temp)

## firstly, map X1 from WD11 to WD19
## (smart_round() ensures total no. of movements match as
## long as grouping is done, even though grouping
## takes a long time)
PW1_new <- PW1 %>%
  left_join(
    select(Ward_Lookup, FID, WD11CD, WD11NM),
    by = c("X1" = "FID")) %>%
  inner_join(
    select(wardweights, wd11cd, wd11nm, wd19cd, wd19nm, areaWeight),
    by = c("WD11CD" = "wd11cd", "WD11NM" = "wd11nm")) %>%
  group_by(X1, X2) %>%
  mutate(X3 = smart_round(X3 * areaWeight)) %>%
  ungroup() %>%
  filter(X3 > 0) %>%
  select(wd19cd, wd19nm, X1, X2, X3)

## check total number of movements out of old X1s match
## after redistributing
stopifnot(group_by(PW1_new, X1) %>% 
  summarise(X3 = sum(X3)) %>%
  inner_join(PlaySize, by = "X1") %>%
  mutate(diff = abs(X3 - X2)) %>%
  pluck("diff") %>%
  {all(. == 0)})

## convert to new X1 based on WD19
PW1_new <- PW1_new %>%
  inner_join(Ward19_Lookup, by = c("wd19cd" = "WD19CD", "wd19nm" = "WD19NM")) %>%
  select(X1 = FID, X2, X3) %>%
  group_by(X1, X2) %>%
  summarise(X3 = sum(X3)) %>%
  ungroup()

## secondly, map X2 from WD11 to WD19
## (smart_round() ensures total no. of movements match as
## long as grouping is done, even though grouping
## takes a long time)
PW1_new <- PW1_new %>%
  left_join(
    select(Ward_Lookup, FID, WD11CD, WD11NM),
    by = c("X2" = "FID")) %>%
  inner_join(
    select(wardweights, wd11cd, wd11nm, wd19cd, wd19nm, areaWeight),
    by = c("WD11CD" = "wd11cd", "WD11NM" = "wd11nm")) %>%
  group_by(X1, X2) %>%
  mutate(X3 = smart_round(X3 * areaWeight)) %>%
  ungroup() %>%
  filter(X3 > 0) %>%
  select(wd19cd, wd19nm, X1, X2, X3)

## check total number of movements into old X2s match
## after redistributing
stopifnot(group_by(PW1_new, X2) %>% 
  summarise(X3 = sum(X3)) %>%
  inner_join(
    PW1 %>%
      group_by(X2) %>%
      summarise(X3 = sum(X3)), 
    by = "X2") %>%
  mutate(diff = abs(X3.x - X3.y)) %>%
  pluck("diff") %>%
  {all(. == 0)})

## convert to new X2 based on WD19
PW1_new <- PW1_new %>%
  inner_join(Ward19_Lookup, by = c("wd19cd" = "WD19CD", "wd19nm" = "WD19NM")) %>%
  select(X1, X2 = FID, X3) %>%
  group_by(X1, X2) %>%
  summarise(X3 = sum(X3)) %>%
  ungroup() %>%
  arrange(X1, X2)

## check total movements match
stopifnot(sum(PW1_new$X3) == sum(PW1$X3))

## now generate new PlaySize
PlaySize_new <- PW1_new %>%
  group_by(X1) %>%
  summarise(X3 = sum(X3)) %>%
  select(X1, X2 = X3) %>%
  arrange(X1)

## create new PlayMatrix as proportions
PlayMatrix_new <- PW1_new %>%
  group_by(X1) %>%
  mutate(X3 = X3 / sum(X3)) %>%
  ungroup() %>%
  arrange(X1, X2)

## create centroid dataset
CBB2019 <- st_centroid(ward19) %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  full_join(Ward19_Lookup, by = c("wd19cd" = "WD19CD", "wd19nm" = "WD19NM"), keep = TRUE)

## check no missing matches
stopifnot(all(!is.na(CBB2019$wd19cd)))
stopifnot(all(!is.na(CBB2019$WD19CD)))

## now order according to the lookup
CBB2019 <- arrange(CBB2019, FID) %>%
  select(X1 = FID, X2 = X, X3 = Y)

## create directory to write out data
pathToNewData <- "2011to2019Data"
dir.create(pathToNewData)

## create description file
description <- c(
"{ \"name\"               : \"2011to2019Data\",",
" \"version\"            : \"Jan 11 2021\",",
" \"author(s)\"          : \"TJ McKinley\",",
"  \"contact(s)\"         : \"t.mckinley@exeter.ac.uk\",",
"  \"work\"               : \"EW19.dat\",",
"  \"work_size\"          : \"WorkSize19.dat\",",
"  \"play\"               : \"PlayMatrix19.dat\",",
"  \"play_size\"          : \"PlaySize19.dat\",",
"  \"position\"           : \"CBB2019.dat\",",
"  \"coordinates\"        : \"x/y\",",
"  \"coordinate_units\"   : \"m\",",
"  \"lookup\"             : \"Ward19_Lookup.csv\",",
"  \"lookup_columns\"     : {\"code\":1, \"name\":2,",
"    \"authority_code\":3, \"authority_name\":4}",
"}")
writeLines(description, paste0(pathToNewData, "/description.json"))

## write all other files
write_delim(EW1_new, paste0(pathToNewData, "/EW19.dat"), col_names = FALSE)
write_delim(WorkSize_new, paste0(pathToNewData, "/WorkSize19.dat"), col_names = FALSE)
write_delim(PlayMatrix_new, paste0(pathToNewData, "/PlayMatrix19.dat"), col_names = FALSE)
write_delim(PlaySize_new, paste0(pathToNewData, "/PlaySize19.dat"), col_names = FALSE)
write_delim(CBB2019, paste0(pathToNewData, "/CBB2019.dat"), col_names = FALSE)
write_csv(Ward19_Lookup, paste0(pathToNewData, "/Ward19_Lookup.csv"))
