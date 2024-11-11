%% INIZIALIZZAZIONE
clear all;
clc;
close all;

%% IMPORTAZIONE DATI
data_circadian_off = load("C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Sim_gabriele_pasti_Circadian_off\Sim_gabriele_pasti_Circadian_off.mat");

%% ESTRAZIONE DATI

%time full simulation
time = data_circadian_off.SimSettings.simDurationMinutes; %[min]

%Gb - Glucose basal rate [mg/dl]
Gb= data_circadian_off.Subjects(1).Params.Gb;

%CR - Ratio carbohydrate-to-insuline ratio [g/U]
CR=data_circadian_off.Subjects(1).Params.CR;

%CF - Correction factor [mg/(dL·U)]
CF=data_circadian_off.Subjects(1).Params.CF;

%Ub_day - daily basal insulin [U/day]
Ub_day = data_circadian_off.Subjects(1).Params.dailyBasalInsulin;

% u(t) - Input insuline [pmol]
inputs_subqInsulin_Normal_Basal_Use= data_circadian_off.Subjects(1).Signals.inputs_subqInsulin_Normal_Basal_Use; 
inputs_subqInsulin_Normal_Bolus_Use=data_circadian_off.Subjects(1).Signals.inputs_subqInsulin_Normal_Bolus_Use; 
u = inputs_subqInsulin_Normal_Basal_Use + inputs_subqInsulin_Normal_Bolus_Use;

%r(t) - rate of carbohydrate(CHO) intake [mg]
r=data_circadian_off.Subjects(1).Signals.inputs_mealCHO; 

%y(t) - sensor [mg/dl]
y = data_circadian_off.Subjects(1).Signals.Sensors__replace_me__;

%% CONVERSIONE UNITA DI MISURA
% r - [mg]->[g]
r = r/1000;
%u(t) - [pmol]->[U]
u = u/6000;
% Ub - [U/day]->[U/min]
Ub = Ub_day/(24*60);


%% STIMA DEI PARAMETRI

%inizializzazione valori iniziali parametri
o1 = 0.004; % [1/min]
o2 = CF; % [mg/(dL·U)]
o3 = CF/CR; % [mg/(dL·g)]
o4 = 56; % [min]
o5 = 40; % [min]
o0 = o1*Gb + o2*Ub; % [mg/dL/min]

%vettore teta0 iniziale
theta0= [o0 o1 o2 o3 o4 o5]';

%divisione dataset train e test
T_4_days=60*24*4;
T_3_days=time-T_4_days;
%definizione tempo di campionamento
Ts = 1;
%stati iniziali
x0=[Gb Ub Ub 0 0]';

min = zeros(6,1);
max = ones(6,1)*inf;
[theta_ott] = fmincon(@(theta) cost_function(theta,x0,y,u,r,Ts,T_4_days,theta0),theta0, [], [], [], [], min, max);


% %% DEFINIZIONE SISTEMA LINEARE
% 
% o0=theta_ott(1);
% o1=theta_ott(2);
% o2=theta_ott(3);
% o3=theta_ott(4);
% o4=theta_ott(5);
% o5=theta_ott(6);
% 
% %Definizione sistema lineare in forma matriciale
% A = [-o1 -o2 0 o3 0;0 -1/o4 1/o4 0 0;0 0 -1/o4 0 0; 0 0 0 -1/o5 1/o5; 0 0 0 0 -1/o5];
% B_u= [0 0 1/o4 0 0]';
% B_r= [0 0 0 0 1/o5]';
% E=[o0 0 0 0 0];
% C=[1 0 0 0 0]';