function print_result(root,data,patient,Ts)

%% COSTRUZIONE PERCORSO
data_root = root + data;

data_real_dir = data_root + "Sim_test_14day_GM\\Sim_test_14day_GM.mat";
data_patients = data_root + "SIM_PAZIENTI\\";

%% ESTRAZIONE DATI REALI

data_real= load(data_real_dir);

[time,Gb,CR,CF,Ub,u,r,y,IOB,Ra] = data_extraction(data_real,patient);

file_id = sprintf('adult#%03d_dati_simulazione_T10081.mat', patient);
file_name = data_patients + file_id;

load(file_name);

o3=theta_ott(4);

% VINCOLO IOB
CHO_UB = 90;
tau = 120;
IOB_s = o3*2*Ub; %[22 - 6)
IOB_d = IOB_s + (CHO_UB/CR + tau*Ub); %[6 - 22)
%tempo di simulazione in minuti
Tmax = 24*60*7;
%costruzione upper-lower bound IOB
IOB_vet = create_IOB_vector(Tmax,IOB_s,IOB_d);

% Indici temporali
n = length(v_x_obs);
t_tot_min = (0:n-1)*Ts; 
startTime = datetime('00:00', 'InputFormat', 'HH:mm');  
t = startTime + minutes(t_tot_min);
t.Format = 'HH:mm'; 

T=(60*24)/Ts;

% salvataggio grafici
output_folder = 'C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\matlab\export_figures\non_linear_control';

%% GRAFICI
figure('Name', ['Risultato simulazione NL - Paziente ' num2str(patient)]);

%---------GLICEMIA---------
ax1= subplot(4, 1, 1);
plot(t(1:end-1),y(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Glicemia reale');
hold on;
plot(t,v_x_obs(1,:), 'b-', 'LineWidth', 1, 'DisplayName', 'Glicemia osservatore (ODO)');
plot(t,G_basal(1:Ts:T*7*Ts+1), 'm--', 'LineWidth', 1, 'DisplayName', 'Gbasal(t)');
xl = xlim;
x_patch = [xl(1) xl(2) xl(2) xl(1)];
y_patch = [70 70 180 180];
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
plot([xl(1) xl(2)], [70 70], 'k--', 'LineWidth', 0.5, 'HandleVisibility','off');
plot([xl(1) xl(2)], [180 180], 'k--', 'LineWidth', 0.5, 'HandleVisibility','off');
uistack(findobj(gca, 'Type', 'line'), 'top');
% Linee verticali per separare i giorni
for day = 1:7
    xline(startTime + days(day), 'k--', 'HandleVisibility', 'off');
end

hold off;
grid on;
xlim([t(1), t(end)]);
ylabel('Glicemia [mg/dL]');
legend('show');
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

%------------IOB------------
IOB_obs = o3*(v_x_obs(2,:) + v_x_obs(3,:));
ax2= subplot(4, 1, 2);
plot(t(1:end-1),IOB_vet(1:Ts:(length(t)-1)*Ts), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
hold on;
plot(t(1:end-1),IOB(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'IOB reale');
plot(t,IOB_obs, 'b-', 'LineWidth', 1, 'DisplayName', 'IOB osservatore (ODO)');
% Linee verticali per separare i giorni
for day = 1:7
    xline(startTime + days(day), 'k--', 'HandleVisibility', 'off');
end
ylabel('IOB [U]');
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

%-----------U in e Pasti-------------
ax3= subplot(4, 1, 3);
yyaxis left;
stem(t,v_u_in(1:Ts:end), 'g', 'filled', 'LineWidth', 1, 'DisplayName', 'Insulina iniettata', 'Marker', 'none');
hold on;
% Linee verticali per separare i giorni
for day = 1:7
    xline(startTime + days(day), 'k--', 'HandleVisibility', 'off');
end
ylabel('Insulina [U/min]');
yyaxis right;
plot(t,rk_in(1:Ts:(length(t))*Ts), 'k-', 'LineWidth', 1, 'DisplayName', 'Pasti');
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
plot(t,Si_tar(1:Ts:(length(t))*Ts), 'm--', 'LineWidth', 1, 'DisplayName', 'S_{I_{Tar}}');
hold on;
plot(t,v_x_obs(6,1:end), 'b-', 'LineWidth', 1, 'DisplayName', 'x_{6_{Obs}}');
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

%Sfondo bianco per l'intera figura
set(gcf, 'Color', 'white','Position', [100, 100, 1200, 800]);

%export fig
% filename = fullfile(output_folder, ['Sim_NL_7days_Patient_' num2str(patient) '.pdf']);
% export_fig(filename, '-pdf');

linkaxes([ax1,ax2,ax3,ax4],'x');

figure('Name', ['EXIT MPC NL - Paziente ' num2str(patient)]);
plot(v_exit);

% disp(['tempo totale: ' num2str(sum(v_time_sol)/60)]);


%% Metriche

disp("PAZIENTE " + num2str(patient));

%media
media = mean(y);              % Calcola la media
percentili = prctile(y, [25 75]); % Calcola il 25° e il 75° percentile

% Mostra i risultati
fprintf('Media Gm: %.2f\n', media);
fprintf('25° percentile: %.2f\n', percentili(1));
fprintf('75° percentile: %.2f\n', percentili(2));

%standard deviation
SD=std(y);
fprintf('SD Gm: %.2f\n', SD);

%coefficienti di variazione CV
CV=SD/media * 100;
fprintf('Coeff. Variazione Gm: %.2f\n', CV);

%tempo tra 70-140
vet_logico_valori_70_140 = (y >= 70) & (y <= 140); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_70_140);
perc_tempo = somma_val_in_range/length(y);
fprintf('Tempo 70-140: %.2f%%\n', perc_tempo * 100);

%tempo tra 70-180
vet_logico_valori_70_180 = (y >= 70) & (y <= 180); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_70_180);
perc_tempo = somma_val_in_range/length(y);
fprintf('Tempo 70-180: %.2f%%\n', perc_tempo * 100);

%tempo tra >180
vet_logico_valori_gt_180 = (y > 180); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_gt_180);
perc_tempo = somma_val_in_range/length(y);
fprintf('Tempo >180: %.2f%%\n', perc_tempo * 100);

%tempo tra >250
vet_logico_valori_gt_250 = (y > 250); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_gt_250);
perc_tempo = somma_val_in_range/length(y);
fprintf('Tempo >250: %.2f%%\n', perc_tempo * 100);

%tempo tra <70
vet_logico_valori_lt_70 = (y < 70); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_lt_70);
num_L1 = somma_val_in_range;
perc_tempo = somma_val_in_range/length(y);
fprintf('Tempo < 70: %.2f%%\n', perc_tempo * 100);

%tempo tra <54
vet_logico_valori_lt_54 = (y < 54); % Vettore logico (1 se nel range, 0 altrimenti)
somma_val_in_range = sum(vet_logico_valori_lt_54);
num_L2 = somma_val_in_range;
perc_tempo = somma_val_in_range/length(y);
fprintf('Tempo < 54: %.2f%%\n', perc_tempo * 100);

%L1 e L2 
fprintf('Num. L1 (<70) : %.2f\n', num_L1);
fprintf('Num. L2 (<54) : %.2f\n', num_L2);

%totale daily insulin
N_day = (24*60);  % Numero di campioni in un giorno
num_days = (length(v_u_in)-1) / N_day;  % Dovrebbe essere 7

tdi_daily_array = zeros(num_days,1);
for i = 1:num_days
    idx_start = (i-1)*N_day + 1;
    idx_end = i*N_day;
    tdi_daily_array(i) = sum(v_u_in(idx_start:idx_end));
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