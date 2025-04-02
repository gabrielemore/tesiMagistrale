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
%modificare numero paziente per selezione
patient = 1;
%parametri iniziali o4 sistema lineare paper Abuin
o4_vet = [56 40 52.5 59.5 45.5 52.5 47.5 50 47.5 50];

%% DIVISIONE TEMPO DI TRAIN E TEST (giorni)
div_train_test_gg_ML = 4; % giorni di train modello lineare
div_train_test_gg_NL = 4; % giorni di train modello non lineare

%% IDENTIFICAZIONE MODELLI
% identificazione modello a 5 stati su dati con circadian OFF
disp('MODELLO 5 STATI - DATASET CIRCADIAN OFF')
theta_ott_ML = model_identification_ML(data_circadian_off_MIO,patient,o4_vet(patient),div_train_test_gg_ML,false,false);
% identificazione modello a 6 stati su dati con circadian ON
disp('MODELLO 6 STATI - DATASET CIRCADIAN ON')
theta_ott_NL = model_identification_NL(data_circadian_on_MIO,patient,div_train_test_gg_NL,theta_ott_ML,true,false);
% identificazione modello a 5 stati su dati con circadian ON
disp('MODELLO 5 STATI - DATASET CIRCADIAN ON')
theta_ott_ML_circ_ON = model_identification_ML(data_circadian_on_MIO,patient,o4_vet(patient),div_train_test_gg_ML,true,false);

%dataset MESSORI
%theta_ott_ML = model_identification_ML(data_circadian_off_MESSORI,patient,o4_vet(patient),div_train_test_gg_ML,true,false);
%theta_ott_NL = model_identification_NL(data_circadian_on_MESSORI,patient,div_train_test_gg_NL,theta_ott_ML,false,false);

%% SALVATAGGIO PARAMETRI

% filename_lineare = sprintf('parametri_LIN_CIRC_ON/theta_ott_lineare_circ_on_p%03d.mat', patient);
% filename_non_lineare = sprintf('parametri_NL/theta_ott_non_lineare_p%03d.mat', patient);
% 
% save(filename_lineare, 'theta_ott_ML_circ_ON');
% save(filename_non_lineare, 'theta_ott_NL');

%% CONFRONTO GRAFICO MODELLO 5 STATI CON MODELLO 6 STATI (con circadian on)
comparison_ML_NL_circON(data_circadian_on_MIO,patient,theta_ott_ML_circ_ON,theta_ott_NL);

%% VALIDAZIONE CON DATASET MESSORI
disp('VALIDAZIONE MESSORI 5 STATI - DATASET CIRCADIAN OFF')
validation_messori_circ_off(data_circadian_off_MESSORI_val,patient,theta_ott_ML,false);
disp('VALIDAZIONE MESSORI 5 STATI - DATASET CIRCADIAN ON')
validation_messori_circ_off(data_circadian_on_MESSORI_val,patient,theta_ott_ML_circ_ON,true);
disp('VALIDAZIONE MESSORI 6 STATI - DATASET CIRCADIAN ON')
validation_messori_circ_on(data_circadian_on_MESSORI_val,patient,theta_ott_NL,false);

%% Metriche identificaizone
metriche_identificazione;

%% VARIBIALI TEMPORALI per il CONTROLLO
Ts = 5;
deltaT = 1;

% %caricamento parametri lineari se precedentemente salvati
% file_id = sprintf('theta_ott_lineare_p%03d.mat', patient);
% file_name = "parametri stimati\\parametri_LIN\\" + file_id;
% load(file_name);
% %caricamento parametri non lineari
% file_id = sprintf('theta_ott_non_lineare_p%03d.mat', patient);
% file_name = "parametri stimati\\parametri_NL\\" + file_id;
% load(file_name);

%% CONTROLLO MODELLO LINEARE
control_linear_model(data_circadian_off_MIO,patient,theta_ott_ML,Ts,deltaT);

%% CONTROLLO MODELLO NON LINEARE
control_NON_linear_model(data_circadian_on_MIO,patient,theta_ott_NL,Ts,deltaT);

%% CONTROLLO MODELLO NON LINEARE v2
control_NON_linear_model_V2(data_circadian_on_MIO,patient,theta_ott_NL,Ts,deltaT);

%% ANALISI SIMULAZIONI UVA/PADOVA
%caricamento dati file simulatore
root_dati_sim = "C:\\Users\\ITAPC\\Documents\\università\\tesi MAGISTRALE\\dati\\";

% ML_PAR_CIRC_OFF = "SIM_14days_PAR_CIRC_OFF\\";
ML_PAR_CIRC_ON = "SIM_14days_PAR_CIRC_ON\\";
NL_PAR_CIRC_ON = "SIM_7days_NL_pesi_rid\\";

%% MODELLO LINEARE con paramentri stimati da modello con CIRCADIAN ON
for i=1:10
    print_result(root_dati_sim,ML_PAR_CIRC_ON,i,Ts);
end
% MEDIA modello lineare con paramentri stimati da modello con CIRCADIAN ON
print_result_mean(root_dati_sim,ML_PAR_CIRC_ON,Ts);
% CVGA modello lineare con paramentri stimati da modello con CIRCADIAN ON
grafici_cvga(root_dati_sim,ML_PAR_CIRC_ON,Ts);

%% MODELLO NON LINEARE
for i=1:10
    print_result_NL(root_dati_sim,NL_PAR_CIRC_ON,i,Ts);
end
% MEDIA modello lineare con paramentri stimati da modello con CIRCADIAN ON
print_result_mean_NL(root_dati_sim,NL_PAR_CIRC_ON,Ts);
% CVGA modello lineare con paramentri stimati da modello con CIRCADIAN ON
grafici_cvga(root_dati_sim,NL_PAR_CIRC_ON,Ts);