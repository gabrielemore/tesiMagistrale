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
validation_messori_circ_off(data_circadian_off_MESSORI_val,patient,theta_ott_ML,true);
disp('VALIDAZIONE MESSORI 5 STATI - DATASET CIRCADIAN ON')
validation_messori_circ_off(data_circadian_on_MESSORI_val,patient,theta_ott_ML_circ_ON,true);
disp('VALIDAZIONE MESSORI 6 STATI - DATASET CIRCADIAN ON')
validation_messori_circ_on(data_circadian_on_MESSORI_val,patient,theta_ott_NL,true);
