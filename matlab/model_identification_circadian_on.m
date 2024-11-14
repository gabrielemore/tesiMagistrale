%% INIZIALIZZAZIONE
clear all;
clc;
close all;

%% IMPORTAZIONE DATI
data_circadian_off = load("C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Sim_gabriele_pasti_Circadian_on\Sim_gabriele_pasti_Circadian_on.mat");

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
y = data_circadian_off.Subjects(1).Signals.Sensors__replace_me__(2:end);

%IOB [U]
IOB_b = data_circadian_off.Subjects(1).Signals.insulin_Normal_Basal_Use_IOB(2:end);
IOB_bol = data_circadian_off.Subjects(1).Signals.insulin_Normal_Bolus_Use_IOB(2:end);
IOB = IOB_b + IOB_bol;

%Rate of glucose apperance in plasma [mg/kg/ min] 
Ra = data_circadian_off.Subjects(1).Signals.subjE_G_app_rate(2:end); 

%% CONVERSIONE UNITA DI MISURA
% r - [mg]->[g]
r = r/1000;
%u(t) - [pmol]->[U]
u = u/6000;
% Ub - [U/day]->[U/min]
Ub = Ub_day/(24*60);

%da verificare
Ra=Ra/2.7;

%% STIMA PARAMETRI INIZIALI 

%inizializzazione valori iniziali parametri
o01 = 0.004; % [1/min]
o02 = CF; % [mg/(dL·U)]
o03 = CF/CR; % [mg/(dL·g)]
o04 = 56; % [min]
o05 = 40; % [min]
o00 = o01*Gb + o02*Ub; % [mg/dL/min]

%vettore teta0 iniziale
theta0= [o00 o01 o02 o03 o04 o05]';

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

%% STIMA PARAMETRI SISTEMA NON LINEARE
%Parametri iniziali per stima parametri sistema non lineare
o0=theta_ott(1);
o1=theta_ott(2);
% o2=theta_ott(3); o2 nel sistema non lineare corrisponde alla sens.
% all'insulina che nel sistema NON lineare non è un parametro fisso
o2=theta_ott(4);
o3=theta_ott(5);
o4=theta_ott(6);

%nuovi parametri inizializzati a 0
o5=1; o6=1; o7=1; o8=1; o9=1;
%vettore teta0 iniziale

theta0= [o0 o1 o2 o3 o4 o5 o6 o7 o8 o9]';

%divisione dataset train e test
T_4_days=60*24*4;
T_3_days=time-T_4_days;
%definizione tempo di campionamento
Ts = 1;
%stati iniziali sistema NON lineare (aggiunto nuovo stato Si)
x0=[Gb Ub Ub 0 0 0]';

min = zeros(10,1);
max = ones(10,1)*inf;
[theta_ott] = fmincon(@(theta) cost_function_non_linear_system(theta,x0,y,u,r,Ts,T_4_days,theta0,CF),theta0, [], [], [], [], min, max);

%%