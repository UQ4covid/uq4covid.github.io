## load libraries
library(tidyverse)
library(lubridate)

## set path
pathToMetaWardsData <- "../../../../../MetaWardsData"
path <- paste0(pathToMetaWardsData, "/model_data/2011Data/")

## load different data sets
CBB2011 <- read_table2(paste0(path, "CBB2011.dat"), col_names = FALSE)
EW1 <- read_table2(paste0(path, "EW1.dat"), col_names = FALSE)
PlayMatrix <- read_table2(paste0(path, "PlayMatrix.dat"), col_names = FALSE)
PlaySize <- read_table2(paste0(path, "PlaySize.dat"), col_names = FALSE)
# seeds <- read_table2(paste0(path, "seeds.dat"), col_names = FALSE) ## not currently used
Ward_Lookup <- read_csv(paste0(path, "Ward_Lookup.csv"), guess_max = 3000)
WorkSize <- read_table2(paste0(path, "WorkSize.dat"), col_names = FALSE)

## CBB2011, Ward_Lookup, WorkSize and PlaySize correspond to a
## ward, indexed by the first column, and some movements of 
## positions in the other columns (I think)

## full details here: https://metawards.org/fileformats/network.html

## check EW1 numbers match to WorkSize
temp <- group_by(EW1, X1) %>%
  summarise(X3 = sum(X3)) %>%
  full_join(WorkSize, by = "X1") %>%
  full_join(PlaySize, by = "X1")
stopifnot(all(temp$X3 == temp$X2.x))
rm(temp)

## check PlayMatrix proportions add up to one across wards
temp <- PlayMatrix %>%
  group_by(X1) %>%
  summarise(X3 = sum(X3)) %>%
  mutate(eq = map_lgl(X3, ~all.equal(., 1)))
stopifnot(all(temp$eq))
rm(temp)

## check for wards with no players
temp_noplay <- EW1 %>%
  select(X1) %>%
  distinct() %>%
  anti_join(
    PlayMatrix %>%
      select(X1) %>%
      distinct(), 
    by = "X1") %>%
  mutate(noplay = 1)

## check generation of PlayMatrix
temp <- full_join(EW1, WorkSize, by = "X1") %>%
  rename(X2 = X2.x) %>%
  mutate(play = X3 / X2.y) %>%
  full_join(PlayMatrix, by = c("X1", "X2")) %>%
  left_join(temp_noplay, by = "X1")
stopifnot(
    filter(temp, !is.na(noplay)) %>%
      pluck("X3.y") %>%
      is.na() %>%
      all()
)
temp <- filter(temp, is.na(noplay))
stopifnot(all.equal(temp$play, temp$X3.y))
rm(temp, temp_noplay)

##################################################
######      MATCH WARDS11 TO WARDS19        ######
##################################################

## check ward lookup
stopifnot(
  group_by(Ward_Lookup, WD11CD) %>%
    count() %>%
    pluck("n") %>%
    {all(. == 1)}
)
stopifnot(all((sort(Ward_Lookup$FID) - 1:nrow(Ward_Lookup)) == 0))

## load Ward19 lookup
Ward19_Lookup <- read_csv("Ward_to_Local_Authority_District_(December_2019)_Lookup_in_the_United_Kingdom.csv")
stopifnot(
  group_by(Ward19_Lookup, WD19CD) %>%
    count() %>%
    pluck("n") %>%
    {all(. == 1)}
)

## read in changes history
changes <- read_csv("Code_History_Database_(Sept_2020)_UK_v2/Changes.csv", guess_max = 300000) %>%
  mutate(OPER_DATE = dmy_hms(OPER_DATE))
changeHist <- read_csv("Code_History_Database_(Sept_2020)_UK_v2/ChangeHistory.csv") %>%
  mutate(OPER_DATE = dmy_hms(OPER_DATE)) %>%
  mutate(TERM_DATE = dmy_hms(TERM_DATE))
equiv <- read_csv("Code_History_Database_(Sept_2020)_UK_v2/Equivalents.csv", guess_max = 300000) %>%
  mutate(OPER_DATE = dmy_hms(OPER_DATE)) %>%
  mutate(TERM_DATE = dmy_hms(TERM_DATE))

## find "live" codes from changesHist and use these to
## extract ENTITYCDs for later filtering
temp <- left_join(
    Ward19_Lookup, 
    filter(changeHist, STATUS == "live"), 
    by = c("WD19CD" = "GEOGCD", "WD19NM" = "GEOGNM"), 
    keep = TRUE)
stopifnot(all(!is.na(temp$GEOGCD)))
stopifnot(all(!is.na(temp$GEOGNM)))
stopifnot(
  group_by(temp, WD19CD, WD19NM) %>%
    count() %>%
    pluck("n") %>%
    {all(. == 1)}
)
## extract ENTITYCDs
entitycds <- unique(temp$ENTITYCD)
changes <- filter(changes, ENTITYCD %in% entitycds)
changeHist <- filter(changeHist, ENTITYCD %in% entitycds)
equiv <- filter(equiv, ENTITYCD %in% entitycds)

## match on ward and LAD name and CD
ward_matches <- inner_join(
    select(Ward_Lookup, WD11NM, WD11CD, LAD11NM, LAD11CD, FID),
    select(Ward19_Lookup, WD19NM, WD19CD, LAD19CD, LAD19NM),
  by = c("WD11CD" = "WD19CD", "WD11NM" = "WD19NM", "LAD11NM" = "LAD19NM", "LAD11CD" = "LAD19CD")) %>%
  mutate(WD19CD = WD11CD, WD19NM = WD11NM, LAD19NM = LAD11NM, LAD19CD = LAD11CD)

## extract wards in 2011 that don't match to 2019
ward_nomatch <- anti_join(
  select(Ward_Lookup, WD11NM, WD11CD, LAD11NM, LAD11CD, FID),
  select(Ward19_Lookup, WD19NM, WD19CD, LAD19CD, LAD19NM),
  by = c("WD11CD" = "WD19CD", "WD11NM" = "WD19NM", "LAD11NM" = "LAD19NM", "LAD11CD" = "LAD19CD"))

## left-join to try to match by WD11CD
temp <- left_join(
    ward_nomatch, 
    changes,
    by = c("WD11CD" = "GEOGCD", "WD11NM" = "GEOGNM"), 
    keep = TRUE) %>%
  arrange(OPER_DATE) %>%
  group_by(WD11CD, WD11NM) %>%
  slice(1) %>%
  ungroup()
## some can't match
filter(temp, is.na(GEOGCD) | is.na(GEOGNM))

## it seems that "Ockbrook And Borrowash" is coded as "Ockbrook & Borrowash"
## in the Ward19 lookup table and "Ockbrook and Borrowash" in the changes table
changes <- mutate(changes, GEOGNM = ifelse(GEOGNM == "Ockbrook And Borrowash", "Ockbrook and Borrowash", GEOGNM)) %>%
  mutate(GEOGNM_P = ifelse(GEOGNM_P == "Ockbrook And Borrowash", "Ockbrook and Borrowash", GEOGNM_P))

## left-join to try to match by WD11CD
## produces a data set where each match is the earliest entry
## which matches WD11CD and WD11NM (missing GEOGCD and GEOGNM rows
## correspond to rows which can't match the "changes" data)
temp <- left_join(
  ward_nomatch, 
  changes, 
  by = c("WD11CD" = "GEOGCD", "WD11NM" = "GEOGNM"), 
  keep = TRUE) %>%
  arrange(OPER_DATE) %>%
  group_by(WD11CD, WD11NM) %>%
  slice(1) %>%
  ungroup()
## check for missing matches
stopifnot(all(!is.na(temp$GEOGCD)))
stopifnot(all(!is.na(temp$GEOGNM)))
## check year range for first entry is <= 2011
stopifnot(all(temp$YEAR <= 2011 | is.na(temp$YEAR)))
## check for single entry per WD11CD/WD11NM combination
stopifnot(
  group_by(temp, WD11CD, WD11NM) %>%
    count() %>%
    pluck("n") %>%
    {all(. == 1)}
)

## now cross-reference WD11CD against "previous" columns in changes table
## and extract most recent change
temp <- left_join(
    select(temp, WD11NM, WD11CD, LAD11NM, LAD11CD, FID),
    changes,
    by = c("WD11CD" = "GEOGCD_P", "WD11NM" = "GEOGNM_P")) %>%
  arrange(OPER_DATE) %>%
  group_by(WD11CD, WD11NM) %>%
  slice(n()) %>%
  ungroup()
## some can't match
filter(temp, is.na(GEOGCD) | is.na(GEOGNM))

## mismatches probably due to not being able to match GEOGCD_P
## and GEOMGNM_P entries, hence check "equiv" database for different 
## codes for these missing values

## extract mismatches
temp1 <- filter(temp, is.na(GEOGCD) | is.na(GEOGNM))

## extract matches (GEOGCD and GEOGNM are now latest ward codes and names)
temp <- filter(temp, !(is.na(GEOGCD) | is.na(GEOGNM)))

## cross-reference mismatches to alternative codes in "equiv" data
## (here GEOGCDO and GEOGNMO are "alternative" codes I think)
temp1 <- left_join(
  select(temp1, WD11NM, WD11CD, LAD11NM, LAD11CD, FID), 
  equiv, 
  by = c("WD11CD" = "GEOGCD", "WD11NM" = "GEOGNM"), 
  keep = TRUE)

## check matches
stopifnot(all(!is.na(temp1$GEOGCDO)))
stopifnot(all(!is.na(temp1$GEOGNMO)))

## now cross-reference "previous" columns in changes table again
## with "equivalent" indexes and extract most recent change
temp1 <- left_join(
  select(temp1, WD11NM, WD11CD, LAD11NM, LAD11CD, FID, GEOGNMO, GEOGCDO),
  changes,
  by = c("GEOGCDO" = "GEOGCD_P", "GEOGNMO" = "GEOGNM_P")) %>%
  arrange(OPER_DATE) %>%
  group_by(WD11CD, WD11NM) %>%
  slice(n()) %>%
  ungroup() %>%
  select(-GEOGCDO, -GEOGNMO)

## some can't match
filter(temp1, is.na(GEOGCD) | is.na(GEOGNM))

## bind matches
temp <- rbind(
  temp, 
  filter(temp1, !(is.na(GEOGCD) | is.na(GEOGNM))) %>%
    select(one_of(colnames(temp)))
)

## extract mismatches
temp1 <- filter(temp1, is.na(GEOGCD) | is.na(GEOGNM))

## mismatches due (I think) to typos in the database
## e.g. "St. Agnes" not "St Agnes" and so forth
filter(changes, GEOGCD_P %in% temp1$WD11CD)
changes <- mutate(changes, GEOGNM = ifelse(GEOGNM == "Looe West and Lansallas", "Looe West and Lansallos", GEOGNM)) %>%
  mutate(GEOGNM_P = ifelse(GEOGNM_P == "Looe West and Lansallas", "Looe West and Lansallos", GEOGNM_P)) %>%
  mutate(GEOGNM = ifelse(GEOGNM == "St Just In Penwith", "St Just in Penwith", GEOGNM)) %>%
  mutate(GEOGNM_P = ifelse(GEOGNM_P == "St Just In Penwith", "St Just in Penwith", GEOGNM_P)) %>%
  mutate(GEOGNM = ifelse(GEOGNM == "St Agnes", "St. Agnes", GEOGNM)) %>%
  mutate(GEOGNM_P = ifelse(GEOGNM_P == "St Agnes", "St. Agnes", GEOGNM_P)) %>%
  mutate(GEOGNM = ifelse(GEOGNM == "St Martin's", "St. Martin's", GEOGNM)) %>%
  mutate(GEOGNM_P = ifelse(GEOGNM_P == "St Martin's", "St. Martin's", GEOGNM_P)) %>%
  mutate(GEOGNM = ifelse(GEOGNM == "St Mary's", "St. Mary's", GEOGNM)) %>%
  mutate(GEOGNM_P = ifelse(GEOGNM_P == "St Mary's", "St. Mary's", GEOGNM_P))

## now try to match again with updated names
temp1 <- left_join(
  select(temp1, WD11NM, WD11CD, LAD11NM, LAD11CD, FID),
  changes,
  by = c("WD11CD" = "GEOGCD_P", "WD11NM" = "GEOGNM_P")) %>%
  arrange(OPER_DATE) %>%
  group_by(WD11CD, WD11NM) %>%
  slice(n()) %>%
  ungroup()

## check matches
stopifnot(all(!is.na(temp1$GEOGCD)))
stopifnot(all(!is.na(temp1$GEOGNM)))

## bind matches
temp <- rbind(temp, temp1)

## cross-reference against Ward19_Lookup
temp <- select(temp, WD11NM, WD11CD, LAD11NM, LAD11CD, GEOGCD, GEOGNM, FID) %>%
  left_join(select(Ward19_Lookup, -FID), by = c("GEOGCD" = "WD19CD", "GEOGNM" = "WD19NM")) %>%
  rename(WD19CD = GEOGCD, WD19NM = GEOGNM) %>%
  select(one_of(colnames(ward_matches)))

## check for missing entries
summarise_all(temp, ~sum(is.na(.)))
filter(temp, is.na(LAD19NM) | is.na(LAD19CD))

## extract matches and mismatches
ward_nomatch <- filter(temp, is.na(LAD19NM) | is.na(LAD19CD))
ward_matches <- rbind(
    ward_matches, 
    filter(temp, !(is.na(LAD19NM) | is.na(LAD19CD)))
  )


GOT TO HERE AND THEN ABANDONED FOR SHAPEFILE APPROACH




## bind to ward_matches
ward_matches <- rbind(ward_matches, temp) %>%
  arrange(FID)
rm(temp, temp1, ward_nomatch)

## check FIDs
stopifnot(all((ward_matches$FID - 1:nrow(ward_matches)) == 0))

## check for no missing values
summarise_all(ward_matches, ~sum(is.na(.)))

filter(ward_matches, is.na(LAD19NM) | is.na(LAD19CD))

## check mismatches manually
as.data.frame(filter(ward_matches, WD11NM != WD19NM & LAD11NM != LAD19NM))

## replace Ward_Lookup
Ward_Lookup <- ward_matches
rm(ward_matches)

##################################################
######        MATCH WARDS TO TIERS          ######
##################################################

## load Tiers data
Tiers_Lookup <- read_csv("tidyLAD19Tiers-2021-01-04.csv") %>%
  select(-codeType)

## extract unique tiers entries
Tiers_Lookup <- arrange(Tiers_Lookup, code, date) %>%
  distinct_at(vars(code, tier), .keep_all = TRUE)

## extract to ward level
Tiers_Lookup <- full_join(Ward_Lookup, Tiers_Lookup, by = c("LAD19CD" = "code", "LAD19NM" = "name"), keep = TRUE)

## check for mismatches
stopifnot(all(!is.na(Tiers_Lookup$LAD19CD)))



GOT TO HERE











## check all wards accounted for at each date
Tiers_Lookup %>%
  group_by(date) %>%
  count()

filter(Tiers_Lookup, is.na(date))


## write file
write.table(CBB2011, "CBB2011.dat", row.names = FALSE, col.names = FALSE)
    
## collapse to LAD level
WorkSize <- left_join(WorkSize, select(Ward_Lookup, ind, LADind), by = c("X1" = "ind"))

## check for missing matches
stopifnot(all(!is.na(WorkSize$LADind)))

## sum in each LAD
WorkSize <- group_by(WorkSize, LADind) %>%
    summarise(X2 = sum(X2)) %>%
    select(LADind, X2)
    
## write file
write.table(WorkSize, "WorkSize.dat", row.names = FALSE, col.names = FALSE)
    
## collapse to LAD level
PlaySize <- left_join(PlaySize, select(Ward_Lookup, ind, LADind), by = c("X1" = "ind"))

## check for missing matches
stopifnot(all(!is.na(PlaySize$LADind)))

## sum in each LAD
PlaySize <- group_by(PlaySize, LADind) %>%
    summarise(X2 = sum(X2)) %>%
    select(LADind, X2)
    
## write file
write.table(PlaySize, "PlaySize.dat", row.names = FALSE, col.names = FALSE)

## I'm figuring EW1 contains movements from ward to ward, with 
## the first two columns denoting wards, and the last denoting
## number of moves

## collapse to LAD level
EW1 <- left_join(EW1, select(Ward_Lookup, ind, LADind), by = c("X1" = "ind")) %>%
    left_join(select(Ward_Lookup, ind, LADind), by = c("X2" = "ind"))

## check for missing matches
stopifnot(all(!is.na(EW1$LADind.x)))
stopifnot(all(!is.na(EW1$LADind.y)))

## sum in each LAD combination
EW1 <- group_by(EW1, LADind.x, LADind.y) %>%
    summarise(X3 = sum(X3)) %>%
    select(LADind.x, LADind.y, X3)
    
## write file
write.table(EW1, "EW1.dat", row.names = FALSE, col.names = FALSE)
  
## PlayMatrix (for some reason) contains proportions instead of counts
## but you can recover EW1 through PlayMatrix and WorkSize
## hence you can recover counts and calculate proportions correctly

PlayMatrix <- left_join(EW1, WorkSize, by = c("LADind.x" = "LADind")) %>%
    left_join(PlaySize, by = c("LADind.x" = "LADind")) %>%
    mutate(X3 = X3 / X2) %>%
    select(-X2)
    
## write file
write.table(PlayMatrix, "PlayMatrix.dat", row.names = FALSE, col.names = FALSE)


