## load libraries
library(tidyverse)
library(patchwork)
library(sf)
library(covid19.nhs.data)

## read in data
deaths <- read_csv("trustDeathsAdmissionsIcuExtract20200624.csv", guess_max = 30000) %>%
    filter(statistic == "death") %>%
    filter(date <= "2020-04-01") %>%
    select(code, date, deaths = value)

## read in data, originally from here: 
## https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=cumDeathsByDeathDate&format=csv
deathsPHE <- read_csv("ltla_2021-05-18.csv", guess_max = 30000) %>%
    filter(date <= "2020-04-01") %>%
    select(code = areaCode, date, cumdeaths = cumDeathsByDeathDate)

## examine cases by date
p1 <- group_by(deathsPHE, date) %>%
    summarise(cumdeaths = sum(cumdeaths), .groups = "drop") %>%
    mutate(data = "PHE") %>%
    rbind(
        group_by(deaths, date) %>%
        summarise(cumdeaths = sum(deaths), .groups = "drop") %>%
        arrange(date) %>%
        mutate(cumdeaths = cumsum(cumdeaths)) %>%
        mutate(data = "linelist")
    ) %>%
    ggplot(aes(x = date, y = cumdeaths, colour = data)) +
        geom_line() + xlab("Date") + ylab("Cumulative Deaths")
ggsave("comp.pdf", p1)

## extract deaths before 14th March 
seeds <- deaths %>%
    arrange(date) %>%
    group_by(date, code) %>%
    summarise(deaths = sum(deaths), .groups = "drop") %>%
    filter(date <= "2020-03-14") %>%
    filter(deaths > 0) %>%
    ungroup()
seedsPHE <- deathsPHE %>%
    arrange(date) %>%
    group_by(date, code) %>%
    summarise(cumdeaths = sum(cumdeaths), .groups = "drop") %>%
    filter(date <= "2020-03-14") %>%
    group_by(code) %>%
    mutate(deaths = cumdeaths - lag(cumdeaths)) %>%
    mutate(deaths = ifelse(is.na(deaths), cumdeaths, deaths)) %>%
    filter(deaths > 0) %>%
    select(-cumdeaths) %>%
    ungroup()

## now read in ward to trust lookup
trust19Lookup <- read_csv("../hospitalCatchments/trust19Lookup.csv")
    
## check trust names match
temp <- anti_join(seeds, trust19Lookup, by = c("code" = "trustId"))

print(paste0(sum(temp$deaths), " DEATHS DON'T MATCH TO TRUSTS IN LOOKUP"))

## temporarily remove mismatches for testing
seeds <- semi_join(seeds, trust19Lookup, by = c("code" = "trustId")) %>%
    mutate(trust = as.numeric(as.factor(code)))

## set path to MetaWardsData
path <- paste0("../../../../../../MetaWardsData/model_data/2011to2019Data/")

## join to Ward19 lookup
Ward19Lookup <- read_csv(paste0(path, "Ward19_Lookup.csv")) %>%
    inner_join(trust19Lookup, by = c("WD19CD" = "code"))

## subset by trusts with initial infections
## using LSHTM lookup
seedsLSHTM <- select(seeds, code, trust) %>%
    distinct() %>%
    inner_join(trust_ltla_mapping, by = c("code" = "trust_code")) %>%
    inner_join(Ward19Lookup, by = c("geo_code" = "LAD19CD")) %>%
    select(WD19CD) %>%
    mutate(seedLSHTM = 1)

## subset by trusts with initial infections
seeds <- select(seeds, code, trust) %>%
    distinct() %>%
    inner_join(Ward19Lookup, by = c("code" = "trustId")) %>%
    select(WD19CD) %>%
    mutate(seed = 1)

## join to Ward19 lookup
seedsPHE <- read_csv(paste0(path, "Ward19_Lookup.csv")) %>%
    inner_join(seedsPHE, by = c("LAD19CD" = "code")) %>%
    select(WD19CD) %>%
    mutate(seedPHE = 1)

## spatial plot of seeds
ward19 <- st_read("../ward11toWard19Mapping/Wards__December_2019__Boundaries_EW_BFC-shp/Wards__December_2019__Boundaries_EW_BFC.shp")

## join seeds to shapefile
ward19 <- left_join(ward19, seeds, by = c("wd19cd" = "WD19CD")) %>%
    left_join(seedsPHE, by = c("wd19cd" = "WD19CD")) %>%
    left_join(seedsLSHTM, by = c("wd19cd" = "WD19CD")) %>%
    filter(!is.na(seed) | !is.na(seedPHE) | !is.na(seedLSHTM))

p1 <- ggplot(filter(ward19, !is.na(seed))) +
    geom_sf() + ggtitle("Line list")
p2 <- ggplot(filter(ward19, !is.na(seedPHE))) +
    geom_sf() + ggtitle("PHE")
p3 <- ggplot(filter(ward19, !is.na(seedLSHTM))) +
    geom_sf() + ggtitle("LSHTM")
p4 <- ggplot(filter(ward19, !is.na(seedPHE) & !is.na(seed) & !is.na(seedLSHTM))) +
    geom_sf() + ggtitle("All")
p <- (p1 + p2) / (p3 + p4)
ggsave("compSpatial.pdf", p, width = 10, height = 7)

st_drop_geometry(ward19) %>%
    inner_join(Ward19Lookup, by = c("wd19cd" = "WD19CD")) %>%
    select(objectid, LAD19NM, seed, seedPHE) %>%
    filter(LAD19NM %in% c("Aylesbury Vale", "Chiltern", "Milton Keynes", "Wokingham", "Wycombe")) %>%
    group_by(LAD19NM, seed, seedPHE) %>%
    count() %>%
    select(-n)

### proportions by region

### load LAD to region lookup table
#regionLookup <- read_csv("Local_Authority_District_to_Region_(April_2019)_Lookup_in_England.csv")

###extract seeds again
#seeds <- deaths %>%
#    arrange(date) %>%
#    group_by(date, code) %>%
#    summarise(deaths = sum(deaths), .groups = "drop") %>%
#    filter(date <= "2020-03-14") %>%
#    filter(deaths > 0) %>%
#    ungroup()
#seedsPHE <- deathsPHE %>%
#    arrange(date) %>%
#    group_by(date, code) %>%
#    summarise(cumdeaths = sum(cumdeaths), .groups = "drop") %>%
#    filter(date <= "2020-03-14") %>%
#    group_by(code) %>%
#    mutate(deaths = cumdeaths - lag(cumdeaths)) %>%
#    mutate(deaths = ifelse(is.na(deaths), cumdeaths, deaths)) %>%
#    filter(deaths > 0) %>%
#    select(-cumdeaths) %>%
#    ungroup()
#    
#seeds <- inner_join(seeds, trust19Lookup, by = c("code" = "trustId"))  %>%
#    inner_join(select(Ward19Lookup, WD19CD, LAD19CD), by = c("code.y" = "WD19CD")) %>%
#    inner_join(select(regionLookup, LAD19CD, RGN19CD, RGN19NM), by = "LAD19CD") %>%
#    group_by(date, RGN19NM) %>%
#    summarise(deaths = sum(deaths))

#p <- deaths %>%
#    inner_join(regionLookup, by = c("areaCode" = "LAD19CD")) %>%
#    group_by(date, RGN19NM) %>%
#    summarise(deaths = sum(cumDeathsByDeathDate), .groups = "drop_last") %>%
#    mutate(prop = deaths / sum(deaths)) %>%
#    ggplot(aes(x = date, y = prop, fill = RGN19NM)) +
#    geom_bar(stat = "identity") + xlab("Date") + 
#    ylab("Proportion cumulative deaths") +
#    labs(fill = "Region")
#ggsave("deathsProps.pdf", p)

