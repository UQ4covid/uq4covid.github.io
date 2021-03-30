## load libraries
library(tidyverse)
library(sf)

## set paths
path <- paste0("../../../../../../MetaWardsData/model_data/2011to2019Data/")

## load catchment shapefile
trust19 <- st_read("WD19toNHSTrust/WD19toNHSTrust.shp")

## remove duplicates
trust19 <- select(trust19, -actBds_, -acutBds) %>%
    distinct(.keep_all = TRUE)

## load point locations of hospitals
trust19Supply <- st_read("WD19toNHSTrustSupply/WD19toNHSTrustSupply.shp")

## load lookup of trusts to wards
trust19Lookup <- read_csv("WD19toNHSTrustCode.csv")

## load Ward19 lookup
Ward19Lookup <- read_csv(paste0(path, "Ward19_Lookup.csv"))

## load Ward19 shapefile
ward19 <- st_read("../Ward11toWard19Mapping/Wards__December_2019__Boundaries_EW_BFC-shp/Wards__December_2019__Boundaries_EW_BFC.shp")

## transform coordinates from WGS84 to British National Grid for trusts
trust19 <- st_transform(trust19, st_crs(ward19))
trust19Supply <- st_transform(trust19Supply, st_crs(ward19))

###########################
######    CHECKS     ######
###########################

## check trusts match with lookup
temp <- full_join(
    st_drop_geometry(trust19) %>%
        select(trustId) %>%
        distinct(),
    st_drop_geometry(trust19Supply) %>%
        select(trustId) %>%
        distinct(),
    keep = TRUE
)
filter(temp, is.na(trustId.x))
filter(temp, is.na(trustId.y))

## remove missing trusts from lookup
trust19Supply <- semi_join(trust19Supply, st_drop_geometry(trust19), by = "trustId")

## plot catchments against locations of hospitals
p <- ggplot() + geom_sf(data = trust19) + geom_sf(data = trust19Supply)
ggsave("catchment_hosp.pdf", p)

#####################################################################
######  MATCH MISSING WARDS VIA CENTROIDS TO CLOSEST HOSPITAL  ######
#####################################################################

## check that all wards accounted for in catchments (should be islands)
mismatches <- anti_join(Ward19Lookup, trust19Lookup, by = c("WD19CD" = "code"))
unique(mismatches$LAD19NM)

## mismatches to wards shapefile
mismatches <- ward19 %>%
    right_join(mismatches, by = c("wd19cd" = "WD19CD", "wd19nm" = "WD19NM"), keep = TRUE)
stopifnot(all(!is.na(mismatches$wd19cd)))

## find minimum distance between wards to match on
## (centroids not used since some trusts split over
## non-contiguous wards with large spatial range)
temp <- st_distance(mismatches, trust19, by_element = FALSE)

## join to nearest hospital and get trust ID
temp <- apply(temp, 1, function(x){
    which(x == min(x))
})
mismatches <- slice(trust19Supply, temp) %>%
    st_drop_geometry() %>%
    select(trustId) %>%
    {cbind(mismatches, .)}

## append mismatches to lookup
trust19Lookup <- rbind(
    trust19Lookup,
    select(st_drop_geometry(mismatches), code = WD19CD, trustId = trustId)
)

## check
mismatches <- anti_join(Ward19Lookup, trust19Lookup, by = c("WD19CD" = "code"))
stopifnot(nrow(mismatches) == 0)

## check that all wards accounted for in catchments (should be islands)
mismatches <- anti_join(trust19Lookup, Ward19Lookup, by = c("code" = "WD19CD"))
stopifnot(nrow(mismatches) == 0)

## save amended trust lookup table
write_csv(trust19Lookup, "trust19Lookup.csv")
write_csv(trust19Lookup, "../../inputs/trust19Lookup.csv")
