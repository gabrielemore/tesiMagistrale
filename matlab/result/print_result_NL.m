function print_result(root,data,patient,Ts)

%% COSTRUZIONE PERCORSO
data_root = root + data;

data_real_dir = data_root + "Sim_test_14day_GM\\Sim_test_14day_GM.mat";
data_patients = data_root + "SIM_PAZIENTI\\";

%% ESTRAZIONE DATI REALI

data_real= load(data_real_dir);

[time,Gb,CR,CF,Ub,u,r,y,IOB,Ra] = data_extraction(data_real,patient);

% file_id = sprintf('adult#%03d_dati_simulazione_T20161.mat', patient);
file_id = sprintf('adult#%03d_dati_simulazione_T1441.mat', patient);
file_name = data_patients + file_id;

load(file_name);

o3=theta_ott(4);

% VINCOLO IOB
CHO_UB = 90;
tau = 120;
IOB_s = o3*2*Ub; %[22 - 6)
IOB_d = IOB_s + (CHO_UB/CR + tau*Ub); %[6 - 22)
%tempo di simulazione in minuti
Tmax = 24*60*1;
%costruzione upper-lower bound IOB
IOB_vet = create_IOB_vector(Tmax,IOB_s,IOB_d);

%% GRAFICI
figure('Name', ['Risultato simulazione NL - Paziente ' num2str(patient)]);

%---------GLICEMIA---------
ax1= subplot(3, 1, 1);
plot(y(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Glicemia reale');
hold on;
plot(v_x_obs(1,:), 'b-', 'LineWidth', 1, 'DisplayName', 'Glicemia osservatore (ODO)');
xl = xlim;
x_patch = [xl(1) xl(2) xl(2) xl(1)];
y_patch = [70 70 140 140];
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
uistack(findobj(gca, 'Type', 'line'), 'top');

hold off;
grid on;
xlim([1, length(v_x_obs)]);
ylabel('Glicemia [mg/dL]');
legend('show');
set(gca, 'FontSize', 12);

%------------IOB------------
IOB_obs = o3*(v_x_obs(2,:) + v_x_obs(3,:));
ax2= subplot(3, 1, 2);
plot(IOB_vet(1:Ts:end), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
hold on;
plot(IOB(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'IOB reale');
plot(IOB_obs, 'b-', 'LineWidth', 1, 'DisplayName', 'IOB osservatore (ODO)');
ylabel('IOB [U]');
hold off;
grid on;
xlim([1, length(IOB_obs)]);
legend('show');
set(gca, 'FontSize', 12);

%-----------U in e Pasti-------------
ax3= subplot(3, 1, 3);
yyaxis left;
stem(v_u_in(1:Ts:end), 'g', 'filled', 'LineWidth', 1, 'DisplayName', 'Insulina iniettata', 'Marker', 'none');
hold on;
ylabel('Insulina [U/min]');
yyaxis right;
plot(rk_in(1:Ts:end), 'k-', 'LineWidth', 1, 'DisplayName', 'Pasti');
ylabel('Ingest [g/min]');
hold off;
grid on;
xlim([0, length(v_x_obs)]);
ylim([0, 5]);
legend('show');
set(gca, 'FontSize', 12);

%Sfondo bianco per l'intera figura
set(gcf, 'Color', 'white');

linkaxes([ax1,ax2,ax3],'x');

figure('Name', ['EXIT MPC NL - Paziente ' num2str(patient)]);
plot(v_exit);

num_var=6;

for i=1:length(v_exit)
    if v_exit(i)==0

        % SOL. PRECEDENTE 
        figure('Name', ['Simulazione PRECEDENTE']);
        e_px1= subplot(2, 1, 1);
        plot(v_xf(num_var*i+1-6,:), 'r-', 'LineWidth', 1, 'DisplayName', 'xf');
        grid on;
        xlim([1, length(v_xf(num_var*i+1-6,:))]);
        legend('show');
        set(gca, 'FontSize', 12);

        e_px2= subplot(2, 1, 2);
        IOB_xf_i = o3*(v_xf(num_var*i+2-6,:) + v_xf(num_var*i+3-6,:));
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
        IOB_xf_i = o3*(v_xf(num_var*i+2,:) + v_xf(num_var*i+3,:));
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
        plot(v_xf(num_var*i+1+6,:), 'r-', 'LineWidth', 1, 'DisplayName', 'xf');
        grid on;
        xlim([1, length(v_xf(num_var*i+1+6,:))]);
        legend('show');
        set(gca, 'FontSize', 12);

        e_pox2= subplot(2, 1, 2);
        IOB_xf_i = o3*(v_xf(num_var*i+2+6,:) + v_xf(num_var*i+3+6,:));
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

