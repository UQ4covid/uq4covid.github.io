## load libraries
library(tidyverse)
library(broom)
library(rgdal)
library(rgeos)
library(maptools)
library(scico)
library(gganimate)
library(gifski)

# ------------------------------------------------------------------------------
# ----------------------------- Read arguments  --------------------------------
# ------------------------------------------------------------------------------

args <- commandArgs(TRUE)
if(length(args) > 0) {
    stopifnot(length(args) == 4)
    ID <- args[1]
    REP <- as.numeric(args[2])
    VAR <- args[3]
    OUTNAME <- args[4]
} else {
  stop("No arguments")
}

## tests
ID <- "Ens0000"
REP <- 1
VAR <- "H"
OUTNAME <- "WeekDeaths"

# ------------------------------------------------------------------------------
# ----------------------- Read shapefiles / lookups ----------------------------
# ------------------------------------------------------------------------------

trust <- readOGR(
  dsn = "shapefiles/WD11-TRUST/WD11-TRUST.shp", 
  stringsAsFactors = FALSE
)

sites <- readOGR(
  dsn = "shapefiles/WD11-TRUST-SITES/WD11-TRUST-SITES.shp", 
  stringsAsFactors = FALSE
) %>% 
  as_tibble() %>%
  rename(long = coords.x1, lat = coords.x2)

trust_df <- broom::tidy(trust, region = "trustId")

ward_lookup <- read_csv("Ward_Lookup.csv") %>%
  dplyr::select(WD11CD, WD11NM, LAD11CD, LAD11NM, ward = FID)

to_trust <- read_csv("WD11ToAcuteTrustIncWalesHB.csv") %>%
  as_tibble()

# ------------------------------------------------------------------------------
# --------------------------------- Read SQL -----------------------------------
# ------------------------------------------------------------------------------

## create ID with replicate attached
ID <- ifelse(REP != 1, paste0(ID, "x", str_pad(REP, 3, pad = "0")), ID)    

## create path
path <- paste0("raw_outputs/", ID, "/stages.db")

## unzip DB
system(paste0("bzip2 -dkf ", path, ".bz2"))

## open connection
con <- DBI::dbConnect(RSQLite::SQLite(), path)
covid <- tbl(con, "compact")

## collect data
data_tmp <- collect(covid) %>%
    dplyr::select(day, ward, !!VAR)
  
## disconnect from DB
DBI::dbDisconnect(con)
  
## remove uncompressed DB
system(paste0("rm ", path))

## join to lookups
plot_tmp <- data_tmp %>%
  rename(var = !!VAR) %>%
  complete(ward = 1:8588, day = 1:max(day), fill = list(var = 0)) %>%
  left_join(ward_lookup, by = "ward") %>%
  left_join(to_trust, by = "WD11CD") %>%
  dplyr::select(
    ward, day, var, trustId
  ) %>%
  group_by(trustId, day) %>%
  summarise(var = mean(var)) %>%
  ungroup() %>%
  right_join(trust_df, by = c("trustId" = "id"))

## set up for gganimate package
p <- ggplot(data = plot_tmp) + 
  geom_polygon( 
    aes(x = long, y = lat, group = group, fill = var)
  ) +
  coord_fixed() +
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_scico(
    paste0("Trust Averaged\n", VAR), 
    direction = 1, 
    palette = "batlow"
  ) +
  transition_time(day) +
  view_follow() +
  labs(title = "Day: {frame_time}")

# Warning: this takes a while...
# If you want a smoother transition mess around with nframes and fps and 
# transformr and ggplot will handle it.
a <- animate(p, nframes = 16, fps = 1, renderer = gifski_renderer())

## write to external files
dir.create("images", showWarnings = FALSE)
anim_save(paste0("images/", OUTNAME, "_", ID, "_", REP, ".gif"), a)

