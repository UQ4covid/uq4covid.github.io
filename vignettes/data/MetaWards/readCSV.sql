-- create tables
create table weeks0(day, week);
create table wards0(ward, week);

-- import lookups
.mode csv weeks0
.import week_lookup.csv weeks0
.mode csv wards0
.import ward_lookup.csv wards0

vacuum;

.exit

