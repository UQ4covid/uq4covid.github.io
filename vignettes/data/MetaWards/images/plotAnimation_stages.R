## load libraries
library(dplyr)
library(readr)
library(tidyr)
library(rgdal)
library(animation)
library(scico)
library(shape)

# ------------------------------------------------------------------------------
# ----------------------------- Read arguments  --------------------------------
# ------------------------------------------------------------------------------

args <- commandArgs(TRUE)
if(length(args) > 0) {
    stopifnot(length(args) == 3)
    ID <- args[1]
    REP <- as.numeric(args[2])
    VAR <- args[3]
} else {
  stop("No arguments")
}

### tests
#ID <- "Ens0000"
#REP <- 1
#VAR <- "H"

print(paste0("ID: ", ID))
print(paste0("REP: ", REP))
print(paste0("VAR: ", VAR))

# ------------------------------------------------------------------------------
# ----------------------- Read shapefiles / lookups ----------------------------
# ------------------------------------------------------------------------------

trust <- readOGR(
    dsn = "../data/WD11-TRUST/WD11-TRUST.shp", 
    stringsAsFactors = FALSE
)

ward_lookup <- read_csv("../data/Ward_Lookup.csv", guess_max = 3000) %>%
    dplyr::select(WD11CD, WD11NM, LAD11CD, LAD11NM, ward = FID)

to_trust <- read_csv("../data/WD11ToAcuteTrustIncWalesHB.csv") %>%
    as_tibble()

# ------------------------------------------------------------------------------
# --------------------------------- Read SQL -----------------------------------
# ------------------------------------------------------------------------------

## create ID with replicate attached
ID <- ifelse(REP != 1, paste0(ID, "x", str_pad(REP, 3, pad = "0")), ID)    

## create path
path <- paste0("../raw_outputs/", ID, "/stages.db")

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
    ungroup() 
    
## set colours for plotting
mycolours <- scico(30, palette = "batlow")
mybreaks <- pretty(seq(0, max(plot_tmp$var), length = 30), n = 30)
plot_tmp <- mutate(plot_tmp, cols = mycolours[findInterval(var, vec = mybreaks)])
mybreaklab <- pretty(mybreaks, n = 5)
mybreaklab <- mybreaklab[mybreaklab <= max(mybreaks)]

## create GIF
saveGIF({
        ## bind subset of data
        for (i in sort(unique(plot_tmp$day))) {
            print(paste("Day", i, "done"))
            ## bind to polygons
            temp <- filter(plot_tmp, day == i) %>%
                dplyr::select(-day)
            temp <- merge(trust, temp, by = "trustId")
            plot(temp, col = temp$cols, main = paste0("Day ", i))
            colorlegend(posx = c(0.82, 0.85), posy = c(0.2, 0.8), 
                zlim = range(mybreaks), zval = mybreaklab, col = mycolours)
        }
    }, 
    movie.name = paste0(ID, "_", REP, "_", VAR, "_day.gif"), 
    ani.width = 500, ani.height = 500, interval = 0.5)
    

