function print_result(data_real,patient,Ts)

%% ESTRAZIONE DATI REALI
[time,Gb,CR,CF,Ub,u,r,y,IOB,Ra] = data_extraction(data_real,patient);

file_id = sprintf('adult#%03d_dati_simulazione_T20161.mat', patient);
file_name = sprintf('C:\\Users\\ITAPC\\Documents\\università\\tesi MAGISTRALE\\dati\\SIM_14days_PAR_CIRC_OFF\\SIM_LIN_PAR_CIRC_OFF\\%s', file_id);

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

%% GRAFICI
figure('Name', ['Risultato simulazione - Paziente ' num2str(patient)]);

%---------GLICEMIA---------
subplot(3, 1, 1);
plot(y(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Glicemia reale');
hold on;
plot(v_x_obs(1,:), 'b-', 'LineWidth', 1, 'DisplayName', 'Glicemia osservatore (ODO)');

% Imposta i limiti dell'asse x per costruire l'area
xl = xlim;

% Crea le coordinate per il patch: l'area compresa tra y=70 e y=140
x_patch = [xl(1) xl(2) xl(2) xl(1)];
y_patch = [70 70 140 140];

% Colora l'area in verde con opacità (trasparenza) 0.3 e nessun bordo
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');

% (Opzionale) Porta il grafico originale in primo piano
uistack(findobj(gca, 'Type', 'line'), 'top');


hold off;
grid on;
xlim([1, length(v_x_obs)]);
ylabel('Glicemia [mg/dL]');
legend('show');
set(gca, 'FontSize', 12);

%------------IOB------------
IOB_obs = o4*(v_x_obs(2,:) + v_x_obs(3,:));
subplot(3, 1, 2);
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
subplot(3, 1, 3);
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

figure('Name', ['EXIT MPC - Paziente ' num2str(patient)]);
plot(v_exit);
end

