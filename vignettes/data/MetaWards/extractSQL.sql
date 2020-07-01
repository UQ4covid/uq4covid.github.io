--.headers on
.mode csv
.import week_lookup.csv weeks
.import ward_lookup.csv wards
create table temp as select compact.*, weeks.week from compact inner join weeks on compact.day = weeks.day;
create table temp1 as select ward, week, sum(H) / 7.0 as Hprev, sum(C) / 7.0 as Cprev, max(DH) + max(DC) as Deaths, :id as output, :rep as replicate from temp group by ward, week;
drop table compact;
drop table temp;
drop table weeks;
create table temp as 
    select * from wards left join temp1 using (ward, week);
drop table temp1;
drop table wards;
update temp set Hprev = 0 where Hprev is null;
update temp set Cprev = 0 where Cprev is null;
update temp set Deaths = 0 where Deaths is null;
update temp set output = :id where output is null;
update temp set replicate = :rep where replicate is null;
alter table temp rename to compact;

