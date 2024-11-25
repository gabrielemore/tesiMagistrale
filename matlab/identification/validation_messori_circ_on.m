%% INIZIALIZZAZIONE
clear all;
clc;
close all;

%% IMPORTAZIONE DATI
%messori
data_circadian_off = load('C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\dati\Messori\ValidationMessori_circ.mat');

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
load('C:\Users\ITAPC\Documents\università\tesi MAGISTRALE\matlab\shared_data\parametri_circadian_on.mat', 'theta_ott_NL');

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

%definizione tempo di campionamento
Ts = 1;
%stati iniziali sistema NON lineare (aggiunto nuovo stato Si)
x0=[Gb Ub Ub 0 0 CF]';

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

%% METRICHE GoF

disp("-------------");
GoF = 100*(1 - (norm(y(1:time)- y_cap(1:time)))/(norm(y(1:time) - mean(y(1:time)))));
disp("GoF messori circ on: " + GoF );
disp("-------------");

%% CONFRONTO GRAFICO

%-------------CGM--------------
figure('Name', ['Confronto CGM: Reale vs Stimato - Paziente ' num2str(patient)]);
% Plot dei dati 
plot(1:1:time, y(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'CGM Reale');
hold on;
plot(1:1:time, y_cap(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'CGM Stimato');

hold off;
grid on;
xlim([0, time]);

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
plot(1:1:time, IOB(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'IOB Reale');
hold on;
plot(1:1:time, IOB_cap(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'IOB Stimato');

hold off;
grid on;
xlim([0, time]);
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
plot(1:1:time, Ra(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Ra Reale');
hold on;
plot(1:1:time, Ra_cap(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Ra Stimato');

hold off;
grid on;
xlim([0, time]);
xlabel('Tempo [min]');
ylabel('Ra [mg/(dl\cdotmin)]');
title(['Confronto Ra: Reale vs Stimato - Paziente ' num2str(patient)]);
legend('show');
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white');