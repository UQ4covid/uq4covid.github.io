## R script to query database
# library(dbplyr)
# library(RSQLite)
library(tidyverse)

## loop over age classes
for(i in 1:8) {
    ## establish connection
    system(paste0("bzip2 -dkf raw_outputs/stages", i, ".db.bz2"))
    con <- DBI::dbConnect(RSQLite::SQLite(), paste0("raw_outputs/stages", i, ".db"))

    ## extract data
    compact <- tbl(con, "compact")

    ## the commands above don't pull down from the database
    ## until you run a 'collect()' pull, but it allows
    ## you to use 'tidyverse'-esque notation, which is converts
    ## to SQL if you're not used to SQL

    ## see e.g. https://db.rstudio.com/dplyr

    ## now extract and reconstruct by filling in gaps
    compact <- collect(compact) %>%
        rename(genpop_0 = Einc, genpop_1 = E, genpop_2 = Iinc, genpop_3 = I, genpop_4 = R, genpop_5 = D) %>%
        rename(asymp_2 = IAinc, asymp_3 = IA, asymp_4 = RA) %>%
        rename(hospital_2 = IHinc, hospital_3 = IH, hospital_4 = RH, hospital_5 = DH) %>%
        rename(critical_2 = ICinc, critical_3 = IC, critical_4 = RC, critical_5 = DC)

    genpop <- select(compact, day, ward, starts_with("genpop_")) %>%
        complete(day = 1:10, nesting(ward)) %>%
        arrange(ward, day) %>%
        mutate_at(vars(starts_with("genpop_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
        fill(starts_with("genpop_"), .direction = "down") %>%
        set_names(c("day", "ward", paste0("stage_", 0:5)))
    asymp <- select(compact, day, ward, starts_with("asymp_")) %>%
        gather(stage, value, -day, -ward) %>%
        complete(stage = paste0("asymp_", 0:5), nesting(ward, day), fill =list(value = 0)) %>%
        spread(stage, value) %>%
        complete(day = 1:10, nesting(ward)) %>%
        arrange(ward, day) %>%
        mutate_at(vars(starts_with("asymp_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
        fill(starts_with("asymp_"), .direction = "down") %>%
        set_names(c("day", "ward", paste0("stage_", 0:5)))
    hospital <- select(compact, day, ward, starts_with("hospital_")) %>%
        gather(stage, value, -day, -ward) %>%
        complete(stage = paste0("hospital_", 0:5), nesting(ward, day), fill =list(value = 0)) %>%
        spread(stage, value) %>%
        complete(day = 1:10, nesting(ward)) %>%
        arrange(ward, day) %>%
        mutate_at(vars(starts_with("hospital_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
        fill(starts_with("hospital_"), .direction = "down") %>%
        set_names(c("day", "ward", paste0("stage_", 0:5)))
    critical <- select(compact, day, ward, starts_with("critical_")) %>%
        gather(stage, value, -day, -ward) %>%
        complete(stage = paste0("critical_", 0:5), nesting(ward, day), fill =list(value = 0)) %>%
        spread(stage, value) %>%
        complete(day = 1:10, nesting(ward)) %>%
        arrange(ward, day) %>%
        mutate_at(vars(starts_with("critical_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
        fill(starts_with("critical_"), .direction = "down") %>%
        set_names(c("day", "ward", paste0("stage_", 0:5)))
    results <- list(GP = genpop, A = asymp, H = hospital, C = critical) %>%
        bind_rows(.id = "demo")

    ## check that there are no individuals in any
    ## class they shouldn't be
    stopifnot(
        filter(results, demo != "GP") %>%
        select(stage_0, stage_1) %>% 
        summarise_all(~all(. == 0)) %>%
        unlist() %>%
        all()
    )

    ## check that there are infections, removals and deaths
    ## except in asymptomatics
    stopifnot(
        filter(results, demo != "A") %>%
        select(demo, stage_2, stage_3, stage_4, stage_5) %>%
        group_by(demo) %>%
        summarise_all(~any(. > 0)) %>%
        gather(stage, value, -demo) %>%
        pluck("value") %>%
        all()
    )
    ## check that there are infections and removals
    ## in asymptomatics
    stopifnot(
        filter(results, demo == "A") %>%
        select(stage_2, stage_3, stage_4) %>%
        summarise_all(~any(. > 0)) %>%
        gather(stage, value) %>%
        pluck("value") %>%
        all()
    )
    ## check for no deaths in asymptomatics
    stopifnot(
        filter(results, demo == "A") %>%
        select(stage_5) %>%
        summarise_all(~all(. == 0)) %>%
        pluck("stage_5")
    )

    ## check for non-decreasing R
    stopifnot(
        group_by(results, demo, ward) %>%
        nest() %>%
        mutate(data = map_lgl(data, function(x) {
            mutate(x, lg = lag(stage_4)) %>%
            mutate(diff = stage_4 - lg) %>%
            pluck("diff") %>% 
            {all(.[!is.na(.)] >=0)}
        })) %>%
        pluck("data") %>%
        all()
    )
    ## check for non-decreasing D
    stopifnot(
        filter(results, demo != "A") %>%
        group_by(demo, ward) %>%
        nest() %>%
        mutate(data = map_lgl(data, function(x) {
            mutate(x, lg = lag(stage_5)) %>%
            mutate(diff = stage_5 - lg) %>%
            pluck("diff") %>% 
            {all(.[!is.na(.)] >=0)}
        })) %>%
        pluck("data") %>%
        all()
    )

    print("genpop")
    filter(genpop, ward == 255)

    print("asymp")
    filter(asymp, ward == 255)

    print("hospital")
    filter(hospital, ward == 255)

    print("critical")
    filter(critical, ward == 255)

    ## remove zero count wards
    ward_counts <- group_by(results, ward) %>%
        summarise(sum = sum(stage_0)) %>%
        filter(sum > 0)
    results <- inner_join(results, ward_counts, by = "ward") %>%
        select(-sum)
        
    ## more comprehensive checks within each individual ward
    ## that has any infection
    temp <- group_by(results, ward) %>%
        nest() %>%
        mutate(data = map_lgl(data, function(x) {
            ## extract incidence
            inc <- select(x, -stage_1, -stage_3)
            colnames(inc)[-c(1:2)] <- c("E", "I", "R", "D")
            inc <- group_by(inc, demo) %>%
                mutate(R = R - lag(R, default = 0), D = D - lag(D, default = 0)) %>%
                ungroup() 
                
            ## extract counts
            counts <- select(x, -stage_0, -stage_2)
            colnames(counts)[-c(1:2)] <- c("E", "I", "R", "D")
            
            ## combine together for tests
            inc <- gather(inc, state, count, -demo, -day) %>%
                filter(demo != "C") %>%
                filter(demo != "H" | state == "I") %>%
                mutate(state = ifelse(demo == "H" & state == "I", "R", state)) %>%
                group_by(day, state) %>%
                summarise_if(is.numeric, sum) %>%
                spread(state, count) %>%
                ungroup() %>%
                mutate(R = R + D) %>%
                select(day, E, I, R)
            
            ## combine together for tests
            counts <- gather(counts, state, count, -demo, -day) %>%
                mutate(state = ifelse(demo == "H" & state != "E", "R", state)) %>%
                mutate(state = ifelse(demo == "C" & state != "E", "R", state)) %>%
                group_by(day, state) %>%
                summarise_if(is.numeric, sum) %>%
                spread(state, count) %>%
                ungroup() %>%
                mutate(R = R + D) %>%
                select(day, E, I, R)
            
            ## check all non-negative
            range_check <- all(c(all(counts >= 0), all(counts >= 0)))
                
            ## check cumulative counts are consistent
            inc_check <- mutate(inc, E = cumsum(E), I = cumsum(I)) %>%
                mutate(d = ifelse(E >= I, 0, 1)) %>%
                mutate(d = d + ifelse(I >= R, 0, 1)) %>%
                pluck("d") %>%
                {all(. == 0)}
                
            ## check we can re-create counts from incidence
            for(i in 2:nrow(inc)) {
                inc$E[i] <- inc$E[i - 1] + inc$E[i] - inc$I[i]
                inc$I[i] <- inc$I[i - 1] + inc$I[i] - inc$R[i]
                inc$R[i] <- inc$R[i - 1] + inc$R[i]
            }
            counts_check <- all(apply(counts - inc, 2, function(x) all(x == 0)))

            ## amalgamate checks
            all(range_check, inc_check, counts_check)
        }))

    stopifnot(all(temp$data))

    print(paste0("All OK - stage ", i))
}
