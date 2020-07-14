
setwd("~/Documents/PostDoc/20200710_covid_plots")

library(dplyr)
library(dbplyr)
library(tidyr)
library(ggplot2)
library(scico)
library(broom)
library(readr)
library(rgdal)
library(DBI)
library(gganimate)
library(transformr)

# ------------------------------------------------------------------------------
# ----------------------------- Read shapefiles --------------------------------
# ------------------------------------------------------------------------------

trust <- readOGR(
  dsn = "shapefiles/WD11-TRUST/WD11-TRUST.shp", 
  stringsAsFactors = F
)

sites <- readOGR(
  dsn = "shapefiles/WD11-TRUST-SITES/WD11-TRUST-SITES.shp", 
  stringsAsFactors = F
) %>% 
  as_tibble() %>%
  rename(long = coords.x1, lat = coords.x2)

trust_df <- broom::tidy(trust, region = "trustId")

# ------------------------------------------------------------------------------
# --------------------------------- Read SQL -----------------------------------
# ------------------------------------------------------------------------------

con <- DBI::dbConnect(RSQLite::SQLite(), 
  "data/uberStages.db"
)

src_dbi(con)

covid <- tbl(con, "compact")

# ------------------------------------------------------------------------------
# --------------------------------- Index SQL ----------------------------------
# - I have already created these indexes locally; you will also have to first --
# ------------------------------------------------------------------------------

# dbSendStatement(con, "CREATE INDEX ward_idx ON compact (ward);")
# dbSendStatement(con, "CREATE INDEX week_idx ON compact (week);")
# dbSendStatement(con, "CREATE INDEX output_idx ON compact (output);")

# ------------------------------------------------------------------------------
# ---------------------------- Plot model output -------------------------------
# ------------------------------------------------------------------------------

ward_lookup <- read_csv("data/ward_lookup.csv") %>%
  dplyr::select(WD11CD, WD11NM, LAD11CD, LAD11NM, ward = FID)

to_trust <- read_csv("data/WD11ToAcuteTrustIncWalesHB.csv") %>%
  as_tibble()

data_tmp <- covid %>%
  filter(output == "Ens0000", replicate == 1) %>%
  collect() %>%
  mutate(
    ward = as.numeric(ward),
    week = as.numeric(week)
  )

plot_tmp <- data_tmp %>%
  left_join(ward_lookup, by = "ward") %>%
  left_join(to_trust, by = "WD11CD") %>%
  dplyr::select(
    ward, week, Hprev, Cprev, Deaths, output, trustId
  ) %>%
  group_by(trustId, week) %>%
  summarise(H_mean = mean(Hprev)) %>%
  ungroup() %>%
  right_join(trust_df, by = c("trustId" = "id"))

p <- ggplot(data = plot_tmp) + 
  geom_polygon( 
    aes(x = long, y = lat, group = group, fill = H_mean)
  ) +
  coord_fixed() +
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_scico(
    "Trust Averaged\nHospital Prevalence", 
    direction = 1, 
    palette = "batlow"
  ) +
  transition_time(week) +
  view_follow() +
  labs(title = "Weeks: {frame_time}")

# Warning: this takes a while...
# If you want a smoother transition mess around with nframes and fps and 
# transformr and ggplot will handle it.
a <- animate(p, nframes = 16, fps = 1) 
anim_save("images/single_sim_batlow.gif", a)

