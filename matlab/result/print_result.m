function print_result(root,data,patient,Ts)

%% COSTRUZIONE PERCORSO
data_root = root + data;

data_real_dir = data_root + "Sim_test_14day_GM\\Sim_test_14day_GM.mat";
data_patients = data_root + "SIM_PAZIENTI\\";

%% ESTRAZIONE DATI REALI

data_real= load(data_real_dir);

[time,Gb,CR,CF,Ub,u,r,y,IOB,Ra] = data_extraction(data_real,patient);

file_id = sprintf('adult#%03d_dati_simulazione_T20161.mat', patient);
file_name = data_patients + file_id;

load(file_name);

% o0=theta_ott(1);
% o1=theta_ott(2);
% o2=theta_ott(3);
% o3=theta_ott(4);
o4=theta_ott(5);
% o5=theta_ott(6);

% VINCOLO IOB
CHO_UB = 90;
tau = 120;
IOB_s = o4*2*Ub; %[22 - 6)
IOB_d = IOB_s + (CHO_UB/CR + tau*Ub); %[6 - 22)
%tempo di simulazione in minuti
Tmax = 24*60*14;
%costruzione upper-lower bound IOB
IOB_vet = create_IOB_vector(Tmax,IOB_s,IOB_d);

% Indici temporali
n = length(v_x_obs);
t_tot_min = (0:n-1)*Ts; 
startTime = datetime('00:00', 'InputFormat', 'HH:mm');  
t = startTime + minutes(t_tot_min);
t.Format = 'HH:mm'; 

%% GRAFICI
figure('Name', ['Risultato simulazione - Paziente ' num2str(patient)]);

%---------GLICEMIA---------
ax1= subplot(3, 1, 1);
plot(t(1:end-1),y(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Glicemia reale');
hold on;
plot(t,v_x_obs(1,:), 'b-', 'LineWidth', 1, 'DisplayName', 'Glicemia osservatore (ODO)');
xl = xlim;
x_patch = [xl(1) xl(2) xl(2) xl(1)];
y_patch = [70 70 180 180];
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
uistack(findobj(gca, 'Type', 'line'), 'top');
% Linee verticali per separare i giorni
for day = 1:14
    xline(startTime + days(day), 'k--', 'HandleVisibility', 'off');
end
hold off;
grid on;
xlim([t(1), t(end)]);
ylabel('Glicemia [mg/dL]');
legend('show');
xtick_vals = [];
for day = 0:13
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
IOB_obs = o4*(v_x_obs(2,:) + v_x_obs(3,:));
ax2= subplot(3, 1, 2);
plot(t(1:end-1),IOB_vet(1:Ts:(length(t)-1)*Ts), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
hold on;
plot(t(1:end-1),IOB(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'IOB reale');
plot(t,IOB_obs, 'b-', 'LineWidth', 1, 'DisplayName', 'IOB osservatore (ODO)');
% Linee verticali per separare i giorni
for day = 1:14
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
ax3= subplot(3, 1, 3);
yyaxis left;
stem(t,v_u_in(1:Ts:end), 'g', 'filled', 'LineWidth', 1, 'DisplayName', 'Insulina iniettata', 'Marker', 'none');
hold on;
% Linee verticali per separare i giorni
for day = 1:14
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

%Sfondo bianco per l'intera figura
set(gcf, 'Color', 'white');

linkaxes([ax1,ax2,ax3],'x');

figure('Name', ['EXIT MPC - Paziente ' num2str(patient)]);
plot(v_exit);

num_var=5;

for i=1:length(v_exit)
    if v_exit(i)==0

        % SOL. PRECEDENTE 
        figure('Name', ['Simulazione PRECEDENTE']);
        e_px1= subplot(2, 1, 1);
        plot(v_xf(num_var*i+1-5,:), 'r-', 'LineWidth', 1, 'DisplayName', 'xf');
        grid on;
        xlim([1, length(v_xf(num_var*i+1-5,:))]);
        legend('show');
        set(gca, 'FontSize', 12);

        e_px2= subplot(2, 1, 2);
        IOB_xf_i = o4*(v_xf(num_var*i+2-5,:) + v_xf(num_var*i+3-5,:));
        plot(IOB_xf_i, 'r-', 'LineWidth', 1, 'DisplayName', 'IOB_xf');
        hold on;
        plot(IOB_vet(i-1:Ts:i+72*Ts), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
        grid on;
        xlim([1, length(IOB_xf_i)]);
        legend('show');
        set(gca, 'FontSize', 12);

        linkaxes([e_px1,e_px2],'x');


        % SOL. NON FATTIBILE
        figure('Name', ['Simulazione NON fattibile']);
        ex1= subplot(2, 1, 1);
        plot(v_xf(num_var*i+1,:), 'r-', 'LineWidth', 1, 'DisplayName', 'xf');
        grid on;
        xlim([1, length(v_xf(num_var*i+1,:))]);
        legend('show');
        set(gca, 'FontSize', 12);

        ex2= subplot(2, 1, 2);
        IOB_xf_i = o4*(v_xf(num_var*i+2,:) + v_xf(num_var*i+3,:));
        plot(IOB_xf_i, 'r-', 'LineWidth', 1, 'DisplayName', 'IOB_xf');
        hold on;
        plot(IOB_vet(i:Ts:i+72*Ts+1), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
        grid on;
        xlim([1, length(IOB_xf_i)]);
        legend('show');
        set(gca, 'FontSize', 12);

        linkaxes([ex1,ex2],'x');

        % SOL. FATTIBILE POST
        figure('Name', ['Simulazione POST']);
        e_pox1= subplot(2, 1, 1);
        plot(v_xf(num_var*i+1+5,:), 'r-', 'LineWidth', 1, 'DisplayName', 'xf');
        grid on;
        xlim([1, length(v_xf(num_var*i+1+5,:))]);
        legend('show');
        set(gca, 'FontSize', 12);

        e_pox2= subplot(2, 1, 2);
        IOB_xf_i = o4*(v_xf(num_var*i+2+5,:) + v_xf(num_var*i+3+5,:));
        plot(IOB_xf_i, 'r-', 'LineWidth', 1, 'DisplayName', 'IOB_xf');
        hold on;
        plot(IOB_vet(i+1:Ts:i+72*Ts+2), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
        grid on;
        xlim([1, length(IOB_xf_i)]);
        legend('show');
        set(gca, 'FontSize', 12);

        linkaxes([e_pox1,e_pox2],'x');

    end

end

%% Metriche

