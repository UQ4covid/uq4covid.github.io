## load libraries
library(tidyverse)
library(magrittr)

## copy data across
system("rm -r test")
system("cp -r /home/tj/Documents/covid/MetaWardsData/model_data/2011to2019Data test")

## read in worksize and playsize
worksize <- read_delim("test/WorkSize19.dat", delim = " ", col_names = FALSE)
playsize <- read_delim("test/PlaySize19.dat", delim = " ", col_names = FALSE)

## aggrgate counts
playsize <- full_join(worksize, playsize, by = "X1") %>%
    mutate(X2 = X2.x + X2.y) %>%
    select(X1, X2) %>%
    arrange(X1)
    
## set zero work counts
worksize$X2 <- 0
EW <- read_delim("test/EW19.dat", delim = " ", col_names = FALSE) %>%
    group_by(X1) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(X3 = 0)
    
# write description
write_delim(worksize, "test/WorkSize19.dat", col_names = FALSE, delim = " ")
write_delim(EW, "test/EW19.dat", col_names = FALSE, delim = " ")

