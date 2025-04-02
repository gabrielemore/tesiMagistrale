function print_result_mean_NL(root,data,Ts)

%% COSTRUZIONE PERCORSO
data_root = root + data;

data_real_dir = data_root + "Sim_test_14day_GM\\Sim_test_14day_GM.mat";
data_patients = data_root + "SIM_PAZIENTI\\";

%% ESTRAZIONE DATI REALI

data_real= load(data_real_dir);

y_mat=[];
IOB_mat=[];
u_in_mat=[];
o3_mat=[];
CR_mat=[];
Ub_mat=[];
Gbasal_mat=[];
Si_mat=[];
x6_obs_mat=[];

for patient=1:10
    if patient ~= 5
    [~,~,CR,~,Ub,~,~,y,IOB,~] = data_extraction(data_real,patient);

    file_id = sprintf('adult#%03d_dati_simulazione_T10081.mat', patient);
    file_name = data_patients + file_id;

    load(file_name);

    % salvo in una matrice
    y_mat= [y_mat; y'];
    IOB_mat = [IOB_mat ;IOB'];
    u_in_mat = [u_in_mat ;v_u_in(1:end-1)'];
    Gbasal_mat = [Gbasal_mat; G_basal];
    Si_mat =[Si_mat;Si_tar];
    x6_obs_mat =[x6_obs_mat;v_x_obs(6,:)];

    %%media per IOB vincolo medio
    o3_mat=[o3_mat theta_ott(4)];
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
mean_U = mean(u_in_mat, 1);         
prct25_U = prctile(u_in_mat, 25, 1);     
prct75_U = prctile(u_in_mat, 75, 1);
%IOB vincolo
o3_mean=mean(o3_mat);
CR_mean=mean(CR_mat);
Ub_mean=mean(Ub_mat);
%Gbasal
Gbasal_mean=mean(Gbasal_mat,1);
%Si_tar
Si_tar_mean=mean(Si_mat,1);
%x6 obs
mean_x6 = mean(x6_obs_mat, 1);         
prct25_x6 = prctile(x6_obs_mat, 25, 1);     
prct75_x6 = prctile(x6_obs_mat, 75, 1);

% VINCOLO IOB
CHO_UB = 90;
tau = 120;
IOB_s = o3_mean*2*Ub_mean; %[22 - 6)
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
output_folder = 'C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\matlab\export_figures\non_linear_control';


%% --- CONFIGURAZIONE FIGURA E VARIABILI DI TEMPO ---
figure('Name', 'Risultato simulazione - Media pazienti');

%% GRAFICO GLICEMIA (Y)
ax1 = subplot(3,1,1);
hold on;

% Plot della media (linea rossa più spessa)
plot(t, mean_Y, 'r-', 'LineWidth', 1, 'DisplayName', 'Media Glicemia');
plot(t,Gbasal_mean(1:length(t)), 'b--', 'LineWidth', 1, 'DisplayName', 'Media Gbasal(t)');
% Patch dell'intervallo target 70-180 (verde, con trasparenza)
xl = xlim;
x_patch = [xl(1) xl(2) xl(2) xl(1)];
y_patch = [70 70 180 180];
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot([xl(1) xl(2)], [70 70], 'k--', 'LineWidth', 0.5, 'HandleVisibility','off');
plot([xl(1) xl(2)], [180 180], 'k--', 'LineWidth', 0.5, 'HandleVisibility','off');

% Aggiungo la banda dei percentili (sfumatura di rosso)
x_coords = [t, fliplr(t)];
y_coords = [prct25_Y, fliplr(prct75_Y)];
% Uso una tonalità di rosso chiaro per l'ombreggiatura
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
stem(t, mean_U, 'g', 'filled', 'LineWidth', 1, 'DisplayName', 'Media U', 'Marker', 'none');

% Aggiungo la banda percentili per U (sfumatura di verde chiaro)
% Uso una tonalità di verde chiaro per la sfumatura: [0.8, 1, 0.8]
x_coords = [t, fliplr(t)];
y_coords = [prct25_U, fliplr(prct75_U)];
fill(x_coords, y_coords, [0.5, 0.9, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility','off');

% Linee verticali per i giorni
for day = 1:7
    xline(startTime + days(day), 'k--', 'HandleVisibility','off');
end
ylabel('Insulina [U/min]');
ylim([0, 4]);
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

%-----------Si_tar e x_obs6-------------
ax4= subplot(4, 1, 4);
plot(t,Si_tar_mean(1:length(t)), 'm--', 'LineWidth', 1, 'DisplayName', 'S_{I_{Tar}}');
hold on;

%interpolazione x6 + percentili
t_obs = linspace(t(1), t(end), length(mean_x6));
x6_obs_interp = interp1(t_obs, mean_x6, t, 'linear');
prct25_x6_interp = interp1(t_obs, prct25_x6, t, 'linear');
prct75_x6_interp = interp1(t_obs, prct75_x6, t, 'linear');


plot(t,x6_obs_interp, 'b-', 'LineWidth', 1, 'DisplayName', 'x_{6_{Obs}}');

x_coords = [t, fliplr(t)];
y_coords = [prct25_x6_interp, fliplr(prct75_x6_interp)];
fill(x_coords, y_coords, [0.5, 0.5, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Linee verticali per separare i giorni
for day = 1:7
    xline(startTime + days(day), 'k--', 'HandleVisibility', 'off');
end
ylabel('S_{I_{Tar}} [mg/(dL·U)]');
hold off;
grid on;
xlim([t(1), t(end)]);
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
linkaxes([ax1,ax2,ax3,ax4],'x');

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
