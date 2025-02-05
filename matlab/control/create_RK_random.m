function [rk,rk_full] = create_RK_random(Tmax)
%NB: al momento Tmax deve essere un multiplo di 1440
%se no,  non viene inserito il pasto nel giorno non completo
MAX_MEAL = 80;
MIN_MEAL = 40;
MAX_SNACK= 20;
MIN_SNACK = 10;

rk=zeros(1,Tmax);
rk_full=zeros(1,Tmax);

%num giorni simulazione
day = floor(Tmax/1440);

for j=1:day
    %COLAZIONE 6-10h (360-600m) (la faccio il 90% delle volte)
    if rand < 0.9
        %inizio un pasto ad un multiplo di 5 (+1 per struttura codice controllo)
        ini=randi([360, 600] / 5) * 5 + 1;
        %durata dle pasto tra 10 e 30 minuti
        duration = randi([10, 30] / 5) * 5;
        fine = ini + duration -1;
        %quantita carbo del pasto
        carbo = randi([MIN_MEAL,MAX_MEAL]);

        %sommo il giorno attuale agli indici inizio e fine
        rk(ini + (j-1)*1440 : fine + (j-1)*1440) = carbo/duration;
        rk_full(ini + (j-1)*1440) = carbo;
    end
    %PRANZO 11.30-14h
    if rand < 0.9
        %inizio un pasto ad un multiplo di 5 (+1 per struttura codice controllo)
        ini=randi([720, 840] / 5) * 5 + 1;
        %durata dle pasto tra 10 e 50 minuti
        duration = randi([10, 50] / 5) * 5;
        fine = ini + duration -1;
        %quantita carbo del pasto
        carbo = randi([MIN_MEAL,MAX_MEAL]);

        %sommo il giorno attuale agli indici inizio e fine
        rk(ini + (j-1)*1440 : fine + (j-1)*1440) = carbo/duration;
        rk_full(ini + (j-1)*1440) = carbo;
    end
    %CENA 19-21h
    if rand < 0.9
        %inizio un pasto ad un multiplo di 5 (+1 per struttura codice controllo)
        ini=randi([1140, 1260] / 5) * 5 + 1;
        %durata dle pasto tra 10 e 50 minuti
        duration = randi([10, 50] / 5) * 5;
        fine = ini + duration -1;
        %quantita carbo del pasto
        carbo = randi([MIN_MEAL,MAX_MEAL]);

        %sommo il giorno attuale agli indici inizio e fine
        rk(ini + (j-1)*1440 : fine + (j-1)*1440) = carbo/duration;
        rk_full(ini + (j-1)*1440) = carbo;
    end
    %SNACK MATTINO 10.45-11.15
    if rand < 0.2
        %inizio un pasto ad un multiplo di 5 (+1 per struttura codice controllo)
        ini=randi([645, 675] / 5) * 5 + 1;
        %durata dle pasto tra 5 e 10 minuti
        duration = randi([5, 10] / 5) * 5;
        fine = ini + duration -1;
        %quantita carbo del pasto
        carbo = randi([MIN_SNACK,MAX_SNACK]);

        %sommo il giorno attuale agli indici inizio e fine
        rk(ini + (j-1)*1440 : fine + (j-1)*1440) = carbo/duration;
        rk_full(ini + (j-1)*1440) = carbo;
    end
    %SNACK POMERIDIANO 16-18
    if rand < 0.4
        %inizio un pasto ad un multiplo di 5 (+1 per struttura codice controllo)
        ini=randi([960, 1080] / 5) * 5 + 1;
        %durata dle pasto tra 10 e 20 minuti
        duration = randi([10, 20] / 5) * 5;
        fine = ini + duration -1;
        %quantita carbo del pasto
        carbo = randi([MIN_SNACK,MAX_SNACK]);

        %sommo il giorno attuale agli indici inizio e fine
        rk(ini + (j-1)*1440 : fine + (j-1)*1440) = carbo/duration;
        rk_full(ini + (j-1)*1440) = carbo;
    end

end
end

