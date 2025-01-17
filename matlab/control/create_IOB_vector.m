function IOB_vet = create_IOB_vector(giorni_in,IOB_s,IOB_d)

% Definizione costanti
MINUTI_GIORNO = 1440;%24 ore
MIN_MATTINA = 360;%6 ore mattina
MIN_GIORNO = 960;%16 ore giorno
MIN_EXTRA = 360;%6 ore extra

%converto in minuti giorni in input
minuti = 24*60*giorni_in;
%Calcolo numero giorni tot
num_giorni = floor(minuti / MINUTI_GIORNO);
%Lunghezza tot vettore IOB 0 minuti + extra 6 ore
minuti_totali = minuti + MIN_EXTRA;
%inizializzo vettore IOB
IOB_vet = zeros(1, minuti_totali);

% Riempimento per ogni giorno completo
for giorno = 0:num_giorni-1
    indice_base = giorno * MINUTI_GIORNO;

    % IOB_s - Mattina (0-6h)
    indice_start = indice_base + 1;
    indice_end = indice_base + MIN_MATTINA;
    IOB_vet(indice_start:indice_end) = IOB_s;

    % IOB_d - Giorno (6-22h)
    indice_start = indice_base + MIN_MATTINA + 1;
    indice_end = indice_base + MIN_MATTINA + MIN_GIORNO;
    IOB_vet(indice_start:indice_end) = IOB_d;

    % Sera (22-24h)
    indice_start = indice_base + MIN_MATTINA + MIN_GIORNO + 1;
    indice_end = indice_base + MINUTI_GIORNO;
    IOB_vet(indice_start:indice_end) = IOB_s;
end

% Extra 6 ore
indice_start = base_idx + MINUTI_GIORNO + 1;
indice_end = base_idx + MINUTI_GIORNO + MIN_EXTRA;
IOB_vet(indice_start:indice_end) = IOB_s;

end

% IOB_vet = zeros(1,1800);
% IOB_vet(1:360) = IOB_s;
% IOB_vet(361:1320) = IOB_d;
% IOB_vet(1321:1440) = IOB_s;
% % + 6h
% IOB_vet(1441:1800) = IOB_s;




% function IOB_vet = create_IOB_vector_by_minutes(total_minutes_input, IOB_s, IOB_d)
%     % Questa funzione crea un vettore IOB per un numero specifico di minuti,
%     % considerando ciclicamente la logica:
%     % - 0h-6h: IOB_s
%     % - 6h-22h: IOB_d
%     % - 22h-24h: IOB_s
%     % E aggiungendo un extra time di 360 minuti alla fine, seguendo la stessa logica.
%     %
%     % Input:
%     % - total_minutes_input: numero totale di minuti (includendo i giorni e/o extra)
%     % - IOB_s: valore da usare per i periodi "s"
%     % - IOB_d: valore da usare per i periodi "d"
%     %
%     % Output:
%     % - IOB_vet: vettore IOB
% 
%     % Durata di un giorno in minuti
%     minutes_per_day = 1440;
% 
%     % Calcolo totale minuti includendo l'extra time di 6 ore
%     total_minutes = total_minutes_input + 360;
% 
%     % Inizializza il vettore IOB
%     IOB_vet = zeros(1, total_minutes);
% 
%     % Riempimento del vettore in base alla logica ciclica
%     for minute = 1:total_minutes
%         % Determina il minuto all'interno del ciclo giornaliero
%         current_day_minute = mod(minute - 1, minutes_per_day) + 1;
% 
%         % Assegna il valore appropriato (IOB_s o IOB_d) in base al minuto
%         if current_day_minute <= 360 || current_day_minute > 1320
%             % Periodo 0h-6h o 22h-24h
%             IOB_vet(minute) = IOB_s;
%         else
%             % Periodo 6h-22h
%             IOB_vet(minute) = IOB_d;
%         end
%     end
% end
