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

### tests
#ID <- "Ens0000"
#REP <- 1
#VAR <- "Deaths"
#OUTNAME <- "WeekDeaths"

# ------------------------------------------------------------------------------
# ----------------------------- Read shapefiles --------------------------------
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

# ------------------------------------------------------------------------------
# --------------------------------- Read SQL -----------------------------------
# ------------------------------------------------------------------------------

con <- DBI::dbConnect(RSQLite::SQLite(), "raw_outputs/uberStages.db")
covid <- tbl(con, "compact")

# ------------------------------------------------------------------------------
# --------------------------------- Index SQL ----------------------------------
# - I have already created these indexes locally; you will also have to first --
# ------------------------------------------------------------------------------

# dbSendStatement(con, "CREATE INDEX output_idx ON compact (output);")

# ------------------------------------------------------------------------------
# ---------------------------- Plot model output -------------------------------
# ------------------------------------------------------------------------------

ward_lookup <- read_csv("Ward_Lookup.csv") %>%
  dplyr::select(WD11CD, WD11NM, LAD11CD, LAD11NM, ward = FID)

to_trust <- read_csv("WD11ToAcuteTrustIncWalesHB.csv") %>%
  as_tibble()

data_tmp <- covid %>%
  filter(output == ID, replicate == REP) %>%
  collect() %>%
  mutate(
    ward = as.numeric(ward),
    week = as.numeric(week)
  )

plot_tmp <- data_tmp %>%
  left_join(ward_lookup, by = "ward") %>%
  left_join(to_trust, by = "WD11CD") %>%
  dplyr::select(
    ward, week, var = !!VAR, output, trustId
  ) %>%
  group_by(trustId, week) %>%
  summarise(var = mean(var)) %>%
  ungroup() %>%
  right_join(trust_df, by = c("trustId" = "id"))

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
  transition_time(week) +
  view_follow() +
  labs(title = "Weeks: {frame_time}")

# Warning: this takes a while...
# If you want a smoother transition mess around with nframes and fps and 
# transformr and ggplot will handle it.
a <- animate(p, nframes = 16, fps = 1, renderer = gifski_renderer())

## write to external files
dir.create("images", showWarnings = FALSE)
anim_save(paste0("images/", OUTNAME, ".gif"), a)

