%% INIZIALIZZAZIONE
clear all;
clc;
close all;

%% IMPORTAZIONE DATI

dataset = 2;

if dataset==1
    %mio
    data_circadian_off = load("C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Sim_gabriele_pasti_Circadian_on\Sim_gabriele_pasti_Circadian_on.mat");
else
    %messori
    data_circadian_off = load('C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Messori\TrainingMessori_circEnabled.mat');
end
%% ESTRAZIONE DATI

patient = 1;

%time full simulation [min]
time = data_circadian_off.SimSettings.simDurationMinutes; 

%Gb - Glucose basal rate [mg/dl]
Gb= data_circadian_off.Subjects(patient).Params.Gb;

%CF - Correction factor [mg/(dL·U)]
CF=data_circadian_off.Subjects(patient).Params.CF;

%Ub_day - daily basal insulin [U/day]
Ub_day = data_circadian_off.Subjects(patient).Params.dailyBasalInsulin;

% u(t) - Input insuline [pmol]
inputs_subqInsulin_Normal_Basal_Use= data_circadian_off.Subjects(patient).Signals.inputs_subqInsulin_Normal_Basal_Use(2:end); 
inputs_subqInsulin_Normal_Bolus_Use=data_circadian_off.Subjects(patient).Signals.inputs_subqInsulin_Normal_Bolus_Use(2:end); 
u = inputs_subqInsulin_Normal_Basal_Use + inputs_subqInsulin_Normal_Bolus_Use;

%r(t) - rate of carbohydrate(CHO) intake [mg]
r=data_circadian_off.Subjects(patient).Signals.inputs_mealCHO(2:end); 

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
%% CARICAMENTO PARAMETRI STIMATI SISTEMA LINEARE
load('parametri_circadian_off.mat', 'theta_ott_ML');

%% STIMA PARAMETRI SISTEMA NON LINEARE
%Parametri iniziali per stima parametri sistema non lineare
o0=theta_ott_ML(1);
o1=theta_ott_ML(2);
% o2=theta_ott(3); o2 nel sistema lineare corrisponde alla sens.
% all'insulina che nel sistema NON lineare non è un parametro fisso
o2=theta_ott_ML(4);
o3=theta_ott_ML(5);
o4=theta_ott_ML(6);

%nuovi parametri
o5=50; o6=25; o7=45; o8=1; o9=1;
%vettore teta0 iniziale

theta0= [o0 o1 o2 o3 o4 o5 o6 o7 o8 o9]';

%divisione dataset train e test
if dataset==1
    %mio
    T_train=60*24*4;
else
    %messori
    T_train = 60*24*3;
end

%definizione tempo di campionamento
Ts = 1;
%stati iniziali sistema NON lineare (aggiunto nuovo stato Si)
x0=[Gb Ub Ub 0 0 CF]';

lb = [0.05,1e-6, 1, 20, 20, 20, 10, 15,  0, 0];
ub = [ 2.3,0.01, 10, 100, 100, 200, 240, 300, 30, 100];

% options = optimoptions('fmincon','Algorithm','interior-point','EnableFeasibilityMode',true, SubproblemAlgorithm='cg');
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxIterations=1000;
options.FunctionTolerance = 1e-5;
options.StepTolerance = 1e-10;
options.MaxFunctionEvaluations=100000;
options.Display='iter-detailed';

[theta_ott_NL] = fmincon(@(theta) cost_function_non_linear_system(theta,x0,y,u,r,Ts,T_train,theta0),theta0, [], [], [], [], lb, ub, [],options);

%% PARAMETRI OTTIMI SISTEMA NON LINEARE

o0=theta_ott_NL(1);
o1=theta_ott_NL(2);
o2=theta_ott_NL(3);
o3=theta_ott_NL(4);
o4=theta_ott_NL(5);
o5=theta_ott_NL(6);
o6=theta_ott_NL(7);
o7=theta_ott_NL(8);
o8=theta_ott_NL(9);
o9=theta_ott_NL(10);

%% SIMULAZIONE SISTEMA NON LINEARE CON PARAMETRI STIMATI

%IOB basal
IOB_basal=2*x0(2)*o3; %x0(2)=Ub

x(:,1) = x0;
y_cap = zeros(time,1);

i = 1:time;
%baseline insulin sensitivity
Si_tar = x0(6)*(1+o8*sin((2*pi*i*Ts)/(60*24) + 2*pi*o9)); % x0(6)=CF
%non abbiamo G_b fisso ma variabile 
G_basal = (o0-Si_tar*x0(2))/o1; %x0(2)=Ub;

for k=1:time
    
    %deviation of gluscose from its basal
    deltaG(k) = x(1,k) - G_basal(k);

    %deviation of IOB from its basal
    IOB(k) = o3*(x(2,k) + x(3,k));
    deltaIOB(k) = IOB(k) - IOB_basal;

    %sistema non lineare
    x_dot(1)= o0 - (o1*x(1,k)) - (x(6,k)*x(2,k)) + (o2*x(4,k)); %G
    x_dot(2)= -(1/o3 * x(2,k)) + (1/o3 * x(3,k)); %Qi
    x_dot(3)= -(1/o3 * x(3,k)) + (1/o3 * u(k)); %Qisub
    x_dot(4)= -(1/o4 * x(4,k)) + (1/o4 * x(5,k)); %Qg
    x_dot(5)= -(1/o4 * x(5,k)) + (1/o4 * r(k)); %Qsto
    x_dot(6)= -(1/o5 * x(6,k)) -(o1/o6 * deltaG(k)) - (1/o7 * deltaIOB(k)) + (1/o5 * Si_tar(k)); %Si

    %discretizzazione
    x(:,k+1) = x(:,k) + Ts*x_dot(:);
    
    y_cap(k) = x(1,k);
   
end

%% METRICHE GoF TRAIN E TEST

disp("-------------");
GoF_train = 100*(1 - (norm(y(1:T_train)- y_cap(1:T_train)))/(norm(y(1:T_train) - mean(y(1:T_train)))));
GoF_test = 100*(1 - (norm(y(T_train+1:time)- y_cap(T_train+1:time)))/(norm(y(T_train+1:time) - mean(y(T_train+1:time)))));
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
IOB_cap = o3*(x(2,:) + x(3,:));

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
Ra_cap = o2*x(4,:);

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
save('parametri_circadian_on.mat', 'theta_ott_NL');
