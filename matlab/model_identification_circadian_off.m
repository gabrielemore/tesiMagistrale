%% INIZIALIZZAZIONE
clear all;
clc;
close all;

%% IMPORTAZIONE DATI
%data_circadian_off = load("C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Sim_gabriele_pasti_Circadian_off\Sim_gabriele_pasti_Circadian_off.mat");

%messori
data_circadian_off = load('C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Messori\TrainingMessori_circDisabled.mat');
%% ESTRAZIONE DATI

patient = 1;

%time full simulation [min]
time = data_circadian_off.SimSettings.simDurationMinutes; 

%Gb - Glucose basal rate [mg/dl]
Gb= data_circadian_off.Subjects(patient).Params.Gb;

%CR - Ratio carbohydrate-to-insuline ratio [g/U]
CR=data_circadian_off.Subjects(patient).Params.CR;

%CF - Correction factor [mg/(dL·U)]
CF=data_circadian_off.Subjects(patient).Params.CF;

%Ub_day - daily basal insulin [U/day]
Ub_day = data_circadian_off.Subjects(patient).Params.dailyBasalInsulin;

% u(t) - Input insuline [pmol]
inputs_subqInsulin_Normal_Basal_Use= data_circadian_off.Subjects(patient).Signals.inputs_subqInsulin_Normal_Basal_Use(2:end); 
inputs_subqInsulin_Normal_Bolus_Use=data_circadian_off.Subjects(patient).Signals.inputs_subqInsulin_Normal_Bolus_Use(2:end); 
u = inputs_subqInsulin_Normal_Basal_Use + inputs_subqInsulin_Normal_Bolus_Use;

%r(t) - rate of carbohydrate(CHO) intake [mg]
r=data_circadian_off.Subjects(11).Signals.inputs_mealCHO(2:end); 

%y(t) - sensor [mg/dl]
y = data_circadian_off.Subjects(patient).Signals.Sensors__replace_me__(2:end);

%IOB [U]
IOB_b = data_circadian_off.Subjects(patient).Signals.insulin_Normal_Basal_Use_IOB(2:end);
IOB_bol = data_circadian_off.Subjects(patient).Signals.insulin_Normal_Bolus_Use_IOB(2:end);
IOB = IOB_b + IOB_bol;

%Rate of glucose apperance in plasma [mg/kg/ min] 
Ra = data_circadian_off.Subjects(patient).Signals.subjE_G_app_rate(2:end); 

%% CONVERSIONE UNITA DI MISURA
% r - [mg]->[g]
r = r/1000;
%u(t) - [pmol]->[U]
u = u/6000;
% Ub - [U/day]->[U/min]
Ub = Ub_day/(24*60);
% Ra
Vg = 2.1143; % [dL/kg];
Ra=Ra/Vg;

%% STIMA DEI PARAMETRI

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
%T_train=60*24*4;
%diviosne dataset Messori
T_train = 60*24*2;
%definizione tempo di campionamento
Ts = 1;
%stati iniziali
x0=[Gb Ub Ub 0 0]';

min = zeros(6,1);
max = ones(6,1)*inf;

% options = optimoptions('fmincon','Algorithm','interior-point','EnableFeasibilityMode',true, SubproblemAlgorithm='cg');
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxIterations=1000;
options.FunctionTolerance = 1e-5;
options.StepTolerance = 1e-10;
options.MaxFunctionEvaluations=100000;
options.Display='iter-detailed';

[theta_ott_ML] = fmincon(@(theta) cost_function(theta,x0,y,u,r,Ts,T_train,theta0),theta0, [], [], [], [], min, max, [],options);

%% DEFINIZIONE SISTEMA LINEARE

o0=theta_ott_ML(1);
o1=theta_ott_ML(2);
o2=theta_ott_ML(3);
o3=theta_ott_ML(4);
o4=theta_ott_ML(5);
o5=theta_ott_ML(6);

%Definizione sistema lineare in forma matriciale
A = [-o1 -o2 0 o3 0;0 -1/o4 1/o4 0 0;0 0 -1/o4 0 0; 0 0 0 -1/o5 1/o5; 0 0 0 0 -1/o5];
B_u= [0 0 1/o4 0 0]';
B_r= [0 0 0 0 1/o5]';
E=[o0 0 0 0 0]';
C=[1 0 0 0 0];

%% CONFRONTO PARAMETRI INIZIALI E PARAMETRI OTTENUTI
CR_cap = o2/o3;
disp("-------------");
disp("CR: " + CR + " - CR stimato: " + CR_cap);
disp("theta4 iniziale: " + o04 + " - theta4 stimato: " + o4);
disp("theta5 iniziale: " + o05 + " - theta5 stimato: " + o5);

%% SIMULAZIONE SISTEMA LINEARE CON PARAMETRI STIMATI

x(:,1)=x0;
for k=1:time
    %sistema lineare
    x(:,k+1) = x(:,k) + Ts*(A*x(:,k) + B_u*u(k) + B_r*r(k) + E);
    y_cap(k)=C*x(:,k);
end

%% METRICHE GoF TRAIN E TEST

disp("-------------");
GoF_train = 100*(1 - (norm(y(1:T_train)- y_cap(1:T_train)'))/(norm(y(1:T_train) - mean(y(1:T_train)))));
GoF_test = 100*(1 - (norm(y(T_train+1:time)- y_cap(T_train+1:time)'))/(norm(y(T_train+1:time) - mean(y(T_train+1:time)))));
disp("GoF train: " + GoF_train );
disp("GoF test: " + GoF_test );

%% CONFRONTO GRAFICO

%-------------CGM--------------
figure('Name', ['Confronto CGM: Reale vs Stimato - Paziente ' num2str(patient)]);
% Plot dei dati 
%train
plot(1:1:T_train, y(1:T_train), 'r--', 'LineWidth', 1.5, 'DisplayName', 'CGM Reale (Train)');
hold on;
plot(1:1:T_train, y_cap(1:T_train), 'b-', 'LineWidth', 1.5, 'DisplayName', 'CGM Stimato (Train)');
%test
plot(T_train+1:1:time, y(T_train+1:time), 'g--', 'LineWidth', 1.5, 'DisplayName', 'CGM Reale (Test)');
plot(T_train+1:1:time, y_cap(T_train+1:time), 'm-', 'LineWidth', 1.5, 'DisplayName', 'CGM Stimato (Test)');

% linea verticale di separazione
xline(T_train, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');

hold off;
grid on;
xlim([0, time]);
% Scritte test e train
ylims = ylim; 
y_pos = ylims(1) + 0.95*(ylims(2)-ylims(1));
text(T_train/2, y_pos, 'Train','HorizontalAlignment', 'center','FontSize', 14,'FontWeight', 'bold');
text(T_train + (time - T_train)/2, y_pos, 'Test','HorizontalAlignment', 'center','FontSize', 14,'FontWeight', 'bold');

xlabel('Tempo [min]');
ylabel('Glucosio [mg/dL]');
title(['Confronto CGM: Reale vs Stimato - Paziente ' num2str(patient)]);
legend('show');
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white');

%------------IOB------------
%calcolo IOB_cap
IOB_cap = o4*(x(2,:) + x(3,:));

figure('Name', ['Confronto IOB: Reale vs Stimato - Paziente ' num2str(patient)]);
% Plot dei dati 
%train
plot(1:1:T_train, IOB(1:T_train), 'r--', 'LineWidth', 1.5, 'DisplayName', 'IOB Reale (Train)');
hold on;
plot(1:1:T_train, IOB_cap(1:T_train), 'b-', 'LineWidth', 1.5, 'DisplayName', 'IOB Stimato (Train)');
%test
plot(T_train+1:1:time, IOB(T_train+1:time), 'g--', 'LineWidth', 1.5, 'DisplayName', 'IOB Reale (Test)');
plot(T_train+1:1:time, IOB_cap(T_train+1:time), 'm-', 'LineWidth', 1.5, 'DisplayName', 'IOB Stimato (Test)');

% linea verticale di separazione
xline(T_train, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');

hold off;
grid on;
xlim([0, time]);
% Scritte test e train
ylims = ylim; 
y_pos = ylims(1) + 0.95*(ylims(2)-ylims(1));
text(T_train/2, y_pos, 'Train','HorizontalAlignment', 'center','FontSize', 14,'FontWeight', 'bold');
text(T_train + (time - T_train)/2, y_pos, 'Test','HorizontalAlignment', 'center','FontSize', 14,'FontWeight', 'bold');

xlabel('Tempo [min]');
ylabel('IOB [U]');
title(['Confronto IOB: Reale vs Stimato - Paziente ' num2str(patient)]);
legend('show');
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white');

%-----------Ra-------------
%calcolo Ra_cap
Ra_cap = o3*x(4,:);

figure('Name', ['Confronto Ra: Reale vs Stimato - Paziente ' num2str(patient)]);
% Plot dei dati 
%train
plot(1:1:T_train, Ra(1:T_train), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Ra Reale (Train)');
hold on;
plot(1:1:T_train, Ra_cap(1:T_train), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Ra Stimato (Train)');
%test
plot(T_train+1:1:time, Ra(T_train+1:time), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ra Reale (Test)');
plot(T_train+1:1:time, Ra_cap(T_train+1:time), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Ra Stimato (Test)');

% linea verticale di separazione
xline(T_train, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');

hold off;
grid on;
xlim([0, time]);
% Scritte test e train
ylims = ylim; 
y_pos = ylims(1) + 0.95*(ylims(2)-ylims(1));
text(T_train/2, y_pos, 'Train','HorizontalAlignment', 'center','FontSize', 14,'FontWeight', 'bold');
text(T_train + (time - T_train)/2, y_pos, 'Test','HorizontalAlignment', 'center','FontSize', 14,'FontWeight', 'bold');

xlabel('Tempo [min]');
ylabel('Ra [mg/(dl\cdotmin)]');
title(['Confronto Ra: Reale vs Stimato - Paziente ' num2str(patient)]);
legend('show');
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white');
%% SALVATAGGIO PARAMETRI SISTEMA LINEARE CIRCADIAN OFF
save('parametri_circadian_off.mat', 'theta_ott_ML');

