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

%% VARIBIALI TEMPORALI per il CONTROLLO
Ts = 5;
deltaT = 1;

%% CONTROLLO MODELLO LINEARE
control_linear_model(data_circadian_off_MIO,patient,theta_ott_ML,Ts,deltaT);

%% CONTROLLO MODELLO NON LINEARE
control_NON_linear_model(data_circadian_on_MIO,patient,theta_ott_NL,Ts,deltaT);

%% GRAFICI
%
% load('dati_intermedi.mat');
% Ts = 5;
% o0=theta_ott_ML(1);
% o1=theta_ott_ML(2);
% o2=theta_ott_ML(3);
% o3=theta_ott_ML(4);
% o4=theta_ott_ML(5);
% o5=theta_ott_ML(6);

% figure('Name', ['Controllo MPC - Paziente ' num2str(patient)]);
%
% %---------GLICEMIA---------
% subplot(5, 1, 1);
% plot(v_y(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Glicemia reale');
% hold on;
% plot(v_xk_obs(1,:), 'b-', 'LineWidth', 1, 'DisplayName', 'Glicemia osservatore (ODO)');
% hold off;
% grid on;
% xlim([1, length(v_xk_obs)]);
% ylabel('Glicemia [mg/dL]');
% legend('show');
% set(gca, 'FontSize', 12);
%
% %------------IOB------------
% IOB_obs = o4*(v_xk_obs(2,:) + v_xk_obs(3,:));
% IOB_real = o4*(v_x_real(2,:) + v_x_real(3,:));
% subplot(5, 1, 2);
% yyaxis left;
% plot(IOB_vet(1:Ts:end), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
% hold on;
% plot(IOB_real(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'IOB reale');
% plot(IOB_obs, 'b-', 'LineWidth', 1, 'DisplayName', 'IOB osservatore (ODO)');
% ylabel('IOB [U]');
% yyaxis right;
% stem(v_u(1:Ts:end), 'g', 'filled', 'LineWidth', 1, 'DisplayName', 'Insulina iniettata', 'Marker', 'none');
% ylabel('Insulina [U/min]');
% hold off;
% grid on;
% xlim([1, length(IOB_obs)]);
% legend('show');
% set(gca, 'FontSize', 12);
%
% %-----------Ra-------------
% Ra_real = o3*v_x_real(4,:);
% Ra_obs = o3*v_xk_obs(4,:);
% subplot(5, 1, 3);
% yyaxis left;
% plot(Ra_real(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Ra reale');
% hold on;
% plot(Ra_obs, 'b-', 'LineWidth', 1, 'DisplayName', 'Ra osservatore (ODO)');
% ylabel('Ra [mg/(dl\cdotmin)]');
% yyaxis right;
% plot(rk_in(1:Ts:end), 'k-', 'LineWidth', 1, 'DisplayName', 'Pasti');
% ylabel('Ingest [g/min]');
% hold off;
% grid on;
% xlim([0, length(v_xk_obs)]);
% ylim([0, 5]);
% legend('show');
% set(gca, 'FontSize', 12);
%
% %-----------D1-------------
% subplot(5, 1, 4);
% plot(v_d1, 'b-', 'LineWidth', 1, 'DisplayName', 'Disturbo');
% grid on;
% xlim([0, length(v_d1)]);
% ylabel('d_1');
% legend('show');
% set(gca, 'FontSize', 12);
%
% %---------COSTO----------
% subplot(5, 1, 5);
% plot(v_VN, 'b-', 'LineWidth', 1, 'DisplayName', 'Costo');
% grid on;
% xlim([0, length(v_VN)]);
% %Etichetta asse X in comune
% xlabel('Istanti [min]');
% ylabel('Vn');
% legend('show');
% set(gca, 'FontSize', 12);
%
% %Sfondo bianco per l'intera figura
% set(gcf, 'Color', 'white');