--.headers on
.mode csv

-- import lookups
.import week_lookup.csv weeks0
.import ward_lookup.csv wards0

-- convert lookups to correct type
alter table weeks0 add column weekI integer;
update weeks0 set weekI = cast(week as integer);
alter table weeks0 add column dayI integer;
update weeks0 set dayI = cast(day as integer);
create table weeks as select dayI as day, weekI as week from weeks0;
drop table weeks0;

alter table wards0 add column weekI integer;
update wards0 set weekI = cast(week as integer);
alter table wards0 add column wardI integer;
update wards0 set wardI = cast(ward as integer);
create table wards as select weekI as week, wardI as ward from wards0;
drop table wards0;

-- join compact to lookup
create table temp as select compact.*, weeks.week from compact inner join weeks on compact.day = weeks.day;
--create table temp1 as select ward, week, sum(H) / 7.0 as Hprev, sum(C) / 7.0 as Cprev, max(DH) + max(DC) as Deaths, :id as output, :rep as replicate from temp group by ward, week;
create table temp1 as select ward, week, sum(H) / 7.0 as Hprev, sum(C) / 7.0 as Cprev, max(DH) + max(DC) as Deaths from temp group by ward, week;
drop table temp;
drop table weeks;
create table temp as 
    select * from wards left join temp1 using (ward, week);
drop table temp1;
drop table wards;
update temp set Hprev = 0 where Hprev is null;
update temp set Cprev = 0 where Cprev is null;
update temp set Deaths = 0 where Deaths is null;
--update temp set output = :id where output is null;
--update temp set replicate = :rep where replicate is null;
alter table temp rename to weeksums;

