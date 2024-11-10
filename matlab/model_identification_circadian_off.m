%% INIZIALIZZAZIONE
clear all;
clc;
close all;

%% IMPORTAZIONE DATI
data_circadian_off = load("C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Sim_gabriele_pasti_Circadian_off\Sim_gabriele_pasti_Circadian_off.mat");

%% ESTRAZIONE DATI

%time full simulation
time = data_circadian_off.SimSettings.simDurationMinutes; %[min]
%patient weight
weigth = data_circadian_off.Subjects(1).Params.BW; %[kg]

%G(t) - Blood glucose concetration [mg/dl]
G=data_circadian_off.Subjects(1).Signals.subjE_Gp_conc; %[mg/dl] Plasma glucose concentration

%Q_i(t) - insulin delivery rate in plasma [U/min] 

%Q_isub(t) - insulin delivery rate in the subcutaneous compartment [U/min]

%u(t) - insulin infusion [U/min]

%Q_g(t) - rate of carbohydrate absorption from the gut [g/min]
Q_g=data_circadian_off.Subjects(1).Signals.subjE_G_util_rate; %[mg/kg/ min] Glucose utilization rate at all tissues 

%Q_sto(t) - glucose delivery rate from the stomach [g/min]
Q_sto=data_circadian_off.Subjects(1).Signals.subjE_G_app_rate; %[mg/kg/ min] Rate of appearance of glucose 

%r(t) - rate of carbohydrate(CHO) intake [g/min]
r=data_circadian_off.Subjects(1).Signals.inputs_mealCHO; %[pmol] meal carbohydrates within the current iteration of the simulation 

% inputs_subqInsulin_Normal_Basal_Use= data_circadian_off.Subjects(1).Signals.inputs_subqInsulin_Normal_Basal_Use; %[pmol] amount of Normal Basal Insulin injected subcutaneously within the current iteration of the simulation   
% inputs_subqInsulin_Normal_Bolus_Use=data_circadian_off.Subjects(1).Signals.inputs_subqInsulin_Normal_Bolus_Use;  %[pmol] amount of Normal Bolus Insulin injected subcutaneously within the current iteration of the simulation  
% u = inputs_subqInsulin_Normal_Basal_Use + inputs_subqInsulin_Normal_Bolus_Use; %[pmol] 

%% CONVERSIONE UNITA DI MISURA

% Q_i [units]->[U/min]
% Q_isub [mg/dl]->[U/min]
% u [pmol]->[U/min]

% Q_g [mg/kg / min]->[g/min]
Q_g = Q_g*weigth/1000;
Q_g=diff(Q_g);
% Q_sto [mg/kg / min]->[g/min]
Q_sto = Q_sto*weigth/1000;
Q_sto = diff(Q_sto);
% r [pmol]->[g/min]
r=r*10^-12*180; %180g/mol massa molare del glucosio
r=diff(r);

%% DEFINIZIONE SISTEMA LINEARE
syms o0 o1 o2 o3 o4 o5; 

A = [-o1 -o2 0 o3 0;0 -1/o4 1/o4 0 0;0 0 -1/o4 0 0; 0 0 0 -1/o5 1/o5; 0 0 0 0 -1/o5];
B_u= [0 0 1/o4 0 0]';
B_r= [0 0 0 0 1/o5]';
E=[o0 0 0 0 0];
C=[1 0 0 0 0]';

%% STIMA DEI PARAMETRI