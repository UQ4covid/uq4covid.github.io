-- convert to correct type
alter table compact add column Einc0 integer;
update compact set Einc0 = cast(Einc as integer);
alter table compact add column E0 integer;
update compact set E0 = cast(E as integer);
alter table compact add column Iinc0 integer;
update compact set Iinc0 = cast(Iinc as integer);
alter table compact add column I0 integer;
update compact set I0 = cast(I as integer);
alter table compact add column RI0 integer;
update compact set RI0 = cast(RI as integer);
alter table compact add column DI0 integer;
update compact set DI0 = cast(DI as integer);
alter table compact add column Ainc0 integer;
update compact set Ainc0 = cast(Ainc as integer);
alter table compact add column A0 integer;
update compact set A0 = cast(A as integer);
alter table compact add column RA0 integer;
update compact set RA0 = cast(RA as integer);
alter table compact add column Hinc0 integer;
update compact set Hinc0 = cast(Hinc as integer);
alter table compact add column H0 integer;
update compact set H0 = cast(H as integer);
alter table compact add column RH0 integer;
update compact set RH0 = cast(RH as integer);
alter table compact add column DH0 integer;
update compact set DH0 = cast(DH as integer);
alter table compact add column Cinc0 integer;
update compact set Cinc0 = cast(Cinc as integer);
alter table compact add column C0 integer;
update compact set C0 = cast(C as integer);
alter table compact add column RC0 integer;
update compact set RC0 = cast(RC as integer);
alter table compact add column DC0 integer;
update compact set DC0 = cast(DC as integer);

create table temp as select day, ward, Einc0 as Einc, E0 as E,
    Iinc0 as Iinc, I0 as I, RI0 as RI, DI0 as DI,
    Ainc0 as Ainc, A0 as A, RA0 as RA,
    Hinc0 as Hinc, H0 as H, RH0 as RH, DH0 as DH,
    Cinc0 as Cinc, C0 as C, RC0 as RC, DC0 as DC from compact;
drop table compact;
alter table temp rename to compact;
vacuum;

