function print_result_mean(root,data,Ts)

%% COSTRUZIONE PERCORSO
data_root = root + data;

data_real_dir = data_root + "Sim_test_14day_GM\\Sim_test_14day_GM.mat";
data_patients = data_root + "SIM_PAZIENTI\\";

%% ESTRAZIONE DATI REALI

data_real= load(data_real_dir);

y_mat=[];
IOB_mat=[];
u_in_mat=[];
o4_mat=[];
CR_mat=[];
Ub_mat=[];

for patient=1:10
    if patient ~= 5
    [~,~,CR,~,Ub,~,~,y,IOB,~] = data_extraction(data_real,patient);

    file_id = sprintf('adult#%03d_dati_simulazione_T20161.mat', patient);
    file_name = data_patients + file_id;

    load(file_name);

    % Considero solo una settimana
    y=y(1:end/2);
    v_x_obs=v_x_obs(:,1:floor(end/2)+1);
    IOB=IOB(1:end/2);
    v_u_in=v_u_in(1:floor(end/2));
    rk_in = rk_in(1:floor(end/2)+1);
    rk_in_full = rk_in_full(1:floor(end/2)+1);
    v_exit =v_exit(1:floor(end/2)+1);

    % salvo in una matrice
    y_mat= [y_mat; y'];
    IOB_mat = [IOB_mat ;IOB'];
    u_in_mat = [u_in_mat ;v_u_in'];

    %%media per IOB vincolo medio
    o4_mat=[o4_mat theta_ott(5)];
    CR_mat=[CR_mat CR];
    Ub_mat=[Ub_mat Ub];
    end
end

%calcolo medie e percentili

%Y
mean_Y = mean(y_mat, 1);         
prct25_Y = prctile(y_mat, 25, 1);     
prct75_Y = prctile(y_mat, 75, 1);
%IOB
mean_IOB = mean(IOB_mat, 1);         
prct25_IOB = prctile(IOB_mat, 25, 1);     
prct75_IOB = prctile(IOB_mat, 75, 1);
%U
mean_U_ML = mean(u_in_mat, 1);         
prct25_U = prctile(u_in_mat, 25, 1);     
prct75_U = prctile(u_in_mat, 75, 1);
%IOB vincolo
o4_mean=mean(o4_mat);
CR_mean=mean(CR_mat);
Ub_mean=mean(Ub_mat);

% VINCOLO IOB
CHO_UB = 90;
tau = 120;
IOB_s = o4_mean*2*Ub_mean; %[22 - 6)
IOB_d = IOB_s + (CHO_UB/CR_mean + tau*Ub_mean); %[6 - 22)
%tempo di simulazione in minuti
Tmax = 24*60*7;
%costruzione upper-lower bound IOB
IOB_vet = create_IOB_vector(Tmax,IOB_s,IOB_d);

% Indici temporali
n = length(mean_Y);
t_tot_min = (0:n-1); 
startTime = datetime('00:00', 'InputFormat', 'HH:mm');  
t = startTime + minutes(t_tot_min);
t.Format = 'HH:mm'; 


% salvataggio grafici
output_folder = 'C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\matlab\export_figures\linear_control';


%% --- CONFIGURAZIONE FIGURA E VARIABILI DI TEMPO ---
figure('Name', 'Risultato simulazione - Media pazienti ');

%% GRAFICO GLICEMIA (Y)
ax1 = subplot(3,1,1);
hold on;

plot(t, mean_Y, 'r-', 'LineWidth', 1, 'DisplayName', 'Media Glicemia');

% Patch dell'intervallo target 70-180 
xl = xlim;
x_patch = [xl(1) xl(2) xl(2) xl(1)];
y_patch = [70 70 180 180];
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot([xl(1) xl(2)], [70 70], 'k--', 'LineWidth', 0.5, 'HandleVisibility','off');
plot([xl(1) xl(2)], [180 180], 'k--', 'LineWidth', 0.5, 'HandleVisibility','off');

% banda dei percentili
x_coords = [t, fliplr(t)];
y_coords = [prct25_Y, fliplr(prct75_Y)];
fill(x_coords, y_coords, [0.9, 0.5, 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility','off');

% Linee verticali per separare i giorni
for day = 1:7
    xline(startTime + days(day), 'k--', 'HandleVisibility','off');
end
hold off;
grid on;
xlim([t(1), t(end)]);
ylim([30, 225]);
ylabel('Glicemia [mg/dL]');
legend('show');

% Impostazioni degli xtick come nel tuo codice
xtick_vals = [];
for day = 0:6
    base = startTime + days(day);
    xtick_vals = [xtick_vals, base + hours(6), base + hours(12), base + hours(18)];
end
xtick_vals = xtick_vals(xtick_vals >= t(1) & xtick_vals <= t(end));
xticks(xtick_vals);
xtickformat('HH:mm');
ax = gca;
tickValues = ax.XTick;
ax.XTickLabel = cellstr(datestr(tickValues, 'HH:MM'));
set(gca, 'FontSize', 12);

%% GRAFICO IOB
ax2 = subplot(3,1,2);
hold on;

% Plot del vincolo IOB (linea nera tratteggiata)
plot(t, IOB_vet(1:length(t)), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
% Plot della media IOB (linea rossa spessa)
plot(t, mean_IOB, 'r-', 'LineWidth', 1, 'DisplayName', 'Media IOB');
% Aggiungo la banda percentili per IOB (sfumatura di rosso chiaro)
x_coords = [t, fliplr(t)];
y_coords = [prct25_IOB, fliplr(prct75_IOB)];
fill(x_coords, y_coords, [0.9, 0.5, 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility','off');

% Linee verticali per i giorni
for day = 1:7
    xline(startTime + days(day), 'k--', 'HandleVisibility','off');
end
ylabel('IOB [U]');
hold off;
grid on;
xlim([t(1), t(end)]);
ylim([0, 13]);
legend('show');
xticks(xtick_vals);
xtickformat('HH:mm');
ax = gca;
tickValues = ax.XTick;
ax.XTickLabel = cellstr(datestr(tickValues, 'HH:MM'));
set(gca, 'FontSize', 12);

%% GRAFICO INSULINA (U)
ax3 = subplot(3,1,3);
hold on;
yyaxis left;
% Plot della media U (linea verde spessa)
stem(t, mean_U, 'g', 'filled', 'LineWidth', 1, 'DisplayName', 'Media u_{in}', 'Marker', 'none');
x_coords = [t, fliplr(t)];
y_coords = [prct25_U, fliplr(prct75_U)];
fill(x_coords, y_coords, [0.5, 0.9, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility','off');
% Linee verticali per i giorni
for day = 1:7
    xline(startTime + days(day), 'k--', 'HandleVisibility','off');
end
ylabel('Insulina [U/min]');
yyaxis right;
plot(t,rk_in(1:length(t)), 'k-', 'LineWidth', 1, 'DisplayName', 'Pasti');
ylabel('Ingest [g/min]');
hold off;
grid on;
xlim([t(1), t(end)]);
ylim([0, 5]);
legend('show');
xticks(xtick_vals);
xtickformat('HH:mm');
ax = gca;
tickValues = ax.XTick;
ax.XTickLabel = cellstr(datestr(tickValues, 'HH:MM'));
set(gca, 'FontSize', 12);

%% Impostazioni finali e salvataggio
set(gcf, 'Color', 'white','Position', [100, 100, 1200, 800]);

% % Esporta in PDF nella cartella desiderata
% filename = fullfile(output_folder, 'Sim_7days_Media.pdf');
% export_fig(filename, '-pdf');

% Link degli assi x tra i subplot
linkaxes([ax1,ax2,ax3],'x');

%% Metriche

disp("PAZIENTE MEDIA");

%media
media = mean(mean_Y);              % Calcola la media
percentili = prctile(mean_Y, [25 75]); % Calcola il 25° e il 75° percentile

% Mostra i risultati
fprintf('Media Gm: %.2f\n', media);
fprintf('25° percentile: %.2f\n', percentili(1));
fprintf('75° percentile: %.2f\n', percentili(2));

%standard deviation
SD=std(mean_Y);
fprintf('SD Gm: %.2f\n', SD);

%coefficienti di variazione CV
CV=SD/media * 100;
fprintf('Coeff. Variazione Gm: %.2f\n', CV);

%tempo tra 70-140
vet_logico_valori_70_140 = (mean_Y >= 70) & (mean_Y <= 140); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_70_140);
perc_tempo = somma_val_in_range/length(mean_Y);
fprintf('Tempo 70-140: %.2f%%\n', perc_tempo * 100);

%tempo tra 70-180
vet_logico_valori_70_180 = (mean_Y >= 70) & (mean_Y <= 180); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_70_180);
perc_tempo = somma_val_in_range/length(mean_Y);
fprintf('Tempo 70-180: %.2f%%\n', perc_tempo * 100);

%tempo tra >180
vet_logico_valori_gt_180 = (mean_Y > 180); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_gt_180);
perc_tempo = somma_val_in_range/length(mean_Y);
fprintf('Tempo >180: %.2f%%\n', perc_tempo * 100);

%tempo tra >250
vet_logico_valori_gt_250 = (mean_Y > 250); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_gt_250);
perc_tempo = somma_val_in_range/length(mean_Y);
fprintf('Tempo >250: %.2f%%\n', perc_tempo * 100);

%tempo tra <70
vet_logico_valori_lt_70 = (mean_Y < 70); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_lt_70);
num_L1 = somma_val_in_range;
perc_tempo = somma_val_in_range/length(mean_Y);
fprintf('Tempo < 70: %.2f%%\n', perc_tempo * 100);

%tempo tra <54
vet_logico_valori_lt_54 = (mean_Y < 54); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_lt_54);
num_L2 = somma_val_in_range;
perc_tempo = somma_val_in_range/length(mean_Y);
fprintf('Tempo < 54: %.2f%%\n', perc_tempo * 100);

%L1 e L2 
fprintf('Num. L1 (<70) : %.2f\n', num_L1);
fprintf('Num. L2 (<54) : %.2f\n', num_L2);

%totale daily insulin
N_day = (24*60);  % Numero di campioni in un giorno
num_days = (length(mean_U)) / N_day;  % Dovrebbe essere 7

tdi_daily_array = zeros(num_days,1);
for i = 1:num_days
    idx_start = (i-1)*N_day + 1;
    idx_end = i*N_day;
    tdi_daily_array(i) = sum(mean_U(idx_start:idx_end));
end

tdi_medio = mean(tdi_daily_array);
percentili = prctile(tdi_daily_array, [25 75]); % Calcola il 25° e il 75° percentile

fprintf('TDI : %.2f\n', tdi_medio);
fprintf('25° percentile: %.2f\n', percentili(1));
fprintf('75° percentile: %.2f\n', percentili(2));


%GMI
GMI = 3.31 + 0.02392 * media;
% Mostra il risultato
fprintf('GMI stimato: %.2f%%\n', GMI);

disp("------------------------------------")
