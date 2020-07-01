--.headers on
attach database :new as newcompact;
insert into main.compact select * from newcompact.compact;
