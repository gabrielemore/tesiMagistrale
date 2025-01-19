%% INIZIALIZZAZIONE
clear;
clc;
close all;

%% SELEZIONE DATASET
%MIO
data_circadian_off_MIO = load("C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Sim_gabriele_pasti_Circadian_off\Sim_gabriele_pasti_Circadian_off.mat");
data_circadian_on_MIO = load("C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Sim_gabriele_pasti_Circadian_on\Sim_gabriele_pasti_Circadian_on.mat");

%MESSORI
data_circadian_off_MESSORI = load('C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Messori\TrainingMessori_circDisabled.mat');
data_circadian_on_MESSORI = load('C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Messori\TrainingMessori_circEnabled.mat');

%MESSORI VALIDATION
data_circadian_off_MESSORI_val = load('C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Messori\ValidationMessori_circDisabled.mat');
data_circadian_on_MESSORI_val = load('C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Messori\ValidationMessori_circ.mat');
%% SELEZIONE PERSONA
patient = 1;
%parametri iniziali o4 sistema lineare paper Abuin
o4_vet = [56 40 52.5 59.5 45.5 52.5 47.5 50 47.5 50];

%% DIVISIONE TEMPO DI TRAIN E TEST (giorni)
div_train_test_gg_ML = 4; % giorni di train modello lineare
div_train_test_gg_NL = 4; % giorni di train modello non lineare

%% IDENTIFICAZIONE MODELLI
%MIO
disp('MODELLO 5 STATI - DATASET CIRCADIAN OFF')
theta_ott_ML = model_identification_ML(data_circadian_off_MIO,patient,o4_vet(patient),div_train_test_gg_ML,false,false);
disp('MODELLO 6 STATI - DATASET CIRCADIAN ON')
theta_ott_NL = model_identification_NL(data_circadian_on_MIO,patient,div_train_test_gg_NL,theta_ott_ML,false,false);
%identificazione modello a 5 stati su dati con circadian ON
disp('MODELLO 5 STATI - DATASET CIRCADIAN ON')
theta_ott_ML_circ_ON = model_identification_ML(data_circadian_on_MIO,patient,o4_vet(patient),div_train_test_gg_ML,false,false);

%dataset MESSORI
%theta_ott_ML = model_identification_ML(data_circadian_off_MESSORI,patient,o4_vet(patient),div_train_test_gg_ML,true,false);
%theta_ott_NL = model_identification_NL(data_circadian_on_MESSORI,patient,div_train_test_gg_NL,theta_ott_ML,false,false);

%% CONFRONTO GRAFICO MODELLO 5 STATI CON MODELLO 6 STATI (con circadian on)
comparison_ML_NL_circON(data_circadian_on_MIO,patient,theta_ott_ML_circ_ON,theta_ott_NL);

%% VALIDAZIONE CON DATASET MESSORI
disp('VALIDAZIONE MESSORI 5 STATI - DATASET CIRCADIAN OFF')
validation_messori_circ_off(data_circadian_off_MESSORI_val,patient,theta_ott_ML,false);
disp('VALIDAZIONE MESSORI 5 STATI - DATASET CIRCADIAN ON')
validation_messori_circ_off(data_circadian_on_MESSORI_val,patient,theta_ott_ML_circ_ON,false);
disp('VALIDAZIONE MESSORI 6 STATI - DATASET CIRCADIAN ON')
validation_messori_circ_on(data_circadian_on_MESSORI_val,patient,theta_ott_NL,false);

%% CONTROLLO MODELLO LINEARE
%T=5 e deltaT=1
control_linear_model(data_circadian_off_MIO,patient,theta_ott_ML,5,1);

% %% GRAFICI
% 
% load('dati_intermedi.mat');
% Ts = 5;
% o0=theta_ott_ML(1);
% o1=theta_ott_ML(2);
% o2=theta_ott_ML(3);
% o3=theta_ott_ML(4);
% o4=theta_ott_ML(5);
% o5=theta_ott_ML(6);
% 
% figure('Name', ['Controllo MPC - Paziente ' num2str(patient)]);
% 
% %---------GLICEMIA---------
% subplot(5, 1, 1); 
% % Plot dei dati 
% plot(v_y(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Glicemia reale');
% hold on;
% plot(v_xk_obs(1,:), 'b-', 'LineWidth', 1, 'DisplayName', 'Glicemia osservatore (ODO)');
% 
% hold off;
% grid on;
% xlim([1, length(v_xk_obs)]);
% ylim([50, 310]);
% xlabel('Istanti [min]');
% ylabel('Glicemia [mg/dL]');
% title(['Glicemia - Paziente ' num2str(patient)]);
% legend('show');
% set(gca, 'FontSize', 12);
% set(gcf, 'Color', 'white');
% 
% %------------IOB------------
% %calcolo IOB_real
% IOB_obs = o4*(v_xk_obs(2,:) + v_xk_obs(3,:));
% IOB_real = o4*(v_x_real(2,:) + v_x_real(3,:));
% 
% % Plot dei dati 
% subplot(5, 1, 2); 
% 
% %asse sinistro
% yyaxis left;
% plot(IOB_vet(1:Ts:end), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
% hold on;
% plot(IOB_real(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'IOB reale');
% plot(IOB_obs, 'b-', 'LineWidth', 1, 'DisplayName', 'IOB osservatore (ODO)');
% ylabel('IOB [U]');
% 
% %asse destro
% yyaxis right;
% stem(v_u(1:Ts:end),'g','filled','LineWidth', 1, 'DisplayName', 'Insulina iniettata','Marker', 'none');
% ylabel('Insulina [U/min]');
% 
% hold off;
% grid on;
% xlim([1, length(IOB_obs)]);
% xlabel('Istanti [min]');
% 
% title(['Confronto IOB: reale vs osservatore vs vincolo - Paziente ' num2str(patient)]);
% legend('show');
% set(gca, 'FontSize', 12);
% set(gcf, 'Color', 'white');
% 
% %-----------Ra-------------
% %calcolo Ra_cap
% Ra_real = o3*v_x_real(4,:);
% Ra_obs = o3*v_xk_obs(4,:);
% 
% % Plot dei dati 
% subplot(5, 1, 3); 
% 
% %asse sinistro
% yyaxis left;
% plot(Ra_real(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Ra reale');
% hold on;
% plot(Ra_obs, 'b-', 'LineWidth', 1, 'DisplayName', 'Ra osservatore (ODO)');
% ylabel('Ra [mg/(dl\cdotmin)]');
% 
% %asse destro
% yyaxis right;
% plot(rk_in(1:Ts:end), 'k-', 'LineWidth', 1, 'DisplayName', 'Pasti');
% ylabel('Ingest [g/min]');
% 
% hold off;
% grid on;
% xlim([0, length(v_xk_obs)]);
% ylim([0, 5]);
% xlabel('Istanti [min]');
% 
% title(['Confronto Ra reale vs osservata - Paziente ' num2str(patient)]);
% legend('show');
% set(gca, 'FontSize', 12);
% set(gcf, 'Color', 'white'); 
% 
% %-----------D1-------------
% % Plot dei dati 
% subplot(5, 1, 4); 
% plot(v_d1, 'b-', 'LineWidth', 1, 'DisplayName', 'Disturbo');
% 
% grid on;
% xlim([0, length(v_d1)]);
% xlabel('Istanti [min]');
% ylabel('d_1');
% title(['Convergenza disturbo d_1 - Paziente ' num2str(patient)]);
% legend('show');
% set(gca, 'FontSize', 12);
% set(gcf, 'Color', 'white'); 
% 
% %---------COSTO----------
% % Plot dei dati 
% subplot(5, 1, 5);
% plot(v_VN, 'b-', 'LineWidth', 1, 'DisplayName', 'Costo');
% 
% grid on;
% xlim([0, length(v_VN)]);
% xlabel('Istanti [min]');
% ylabel('Vn');
% title(['Costo - Paziente ' num2str(patient)]);
% legend('show');
% set(gca, 'FontSize', 12);
% set(gcf, 'Color', 'white'); 
