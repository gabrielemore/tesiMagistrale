function comparison_ML_NL_circON(data_circadian_on,patient,theta_ott_ML,theta_ott_NL)
%% ESTRAZIONE DATI
[time,Gb,~,CF,Ub,u,r,y,IOB,Ra] = data_extraction(data_circadian_on,patient);

%% INIZIALIZZAZIONI SIMULAZIONE
%definizione tempo di campionamento
Ts = 1;
%stati iniziali
x0_ML=[Gb Ub Ub 0 0]';
%stati iniziali sistema NON lineare (aggiunto nuovo stato Si)
x0_NL=[Gb Ub Ub 0 0 CF]';

%% DEFINIZIONE SISTEMA LINEARE
o0_ML=theta_ott_ML(1);
o1_ML=theta_ott_ML(2);
o2_ML=theta_ott_ML(3);
o3_ML=theta_ott_ML(4);
o4_ML=theta_ott_ML(5);
o5_ML=theta_ott_ML(6);

%Definizione sistema lineare in forma matriciale
A = [-o1_ML -o2_ML 0 o3_ML 0;0 -1/o4_ML 1/o4_ML 0 0;0 0 -1/o4_ML 0 0; 0 0 0 -1/o5_ML 1/o5_ML; 0 0 0 0 -1/o5_ML];
B_u= [0 0 1/o4_ML 0 0]';
B_r= [0 0 0 0 1/o5_ML]';
E=[o0_ML 0 0 0 0]';
C=[1 0 0 0 0];

%% SIMULAZIONE SISTEMA LINEARE CON PARAMETRI STIMATI
y_cap_ML = zeros(time,1);
x_ML(:,1)=x0_ML;
for k=1:time
    %sistema lineare
    x_ML(:,k+1) = x_ML(:,k) + Ts*(A*x_ML(:,k) + B_u*u(k) + B_r*r(k) + E);
    y_cap_ML(k)=C*x_ML(:,k);
end

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
IOB_basal=2*x0_NL(2)*o3; %x0_NL(2)=Ub

x_NL(:,1) = x0_NL;
y_cap_NL = zeros(time,1);
deltaG = zeros(time,1);
deltaIOB = zeros(time,1);
IOB_NL = zeros(time,1);

i = 1:time;
%baseline insulin sensitivity
Si_tar = x0_NL(6)*(1+0.01*o8*sin((2*pi*i*Ts)/(60*24) + 2*pi*0.01*o9)); % x0_NL(6)=CF
%non abbiamo G_b fisso ma variabile
G_basal = (o0-Si_tar*x0_NL(2))/o1; %x0(2)=Ub;

for k=1:time
    %deviation of gluscose from its basal
    deltaG(k) = x_NL(1,k) - G_basal(k);

    %deviation of IOB from its basal
    IOB_NL(k) = o3*(x_NL(2,k) + x_NL(3,k));
    deltaIOB(k) = IOB_NL(k) - IOB_basal;

    %sistema non lineare
    x_dot(1)= o0 - (o1*x_NL(1,k)) - (x_NL(6,k)*x_NL(2,k)) + (o2*x_NL(4,k)); %G
    x_dot(2)= -(1/o3 * x_NL(2,k)) + (1/o3 * x_NL(3,k)); %Qi
    x_dot(3)= -(1/o3 * x_NL(3,k)) + (1/o3 * u(k)); %Qisub
    x_dot(4)= -(1/o4 * x_NL(4,k)) + (1/o4 * x_NL(5,k)); %Qg
    x_dot(5)= -(1/o4 * x_NL(5,k)) + (1/o4 * r(k)); %Qsto
    x_dot(6)= -(1/o5 * x_NL(6,k)) -(o1/o6 * deltaG(k)) - (1/o7 * deltaIOB(k)) + (1/o5 * Si_tar(k)); %Si

    %discretizzazione
    x_NL(:,k+1) = x_NL(:,k) + Ts*x_dot(:);

    y_cap_NL(k) = x_NL(1,k);

end

%% CONFRONTO GRAFICO
%-------------CGM--------------
figure('Name', ['Confronto ML vs NL lineare - Paziente ' num2str(patient)]);
subplot(3, 1, 1); 
% Plot dei dati 
plot(1:1:time, y(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'CGM reale');
hold on;
plot(1:1:time, y_cap_NL(1:time), 'k-', 'LineWidth', 1.5, 'DisplayName', 'CGM NL');
plot(1:1:time, y_cap_ML(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'CGM ML');

hold off;
grid on;
xlim([0, time]);

xlabel('Tempo [min]');
ylabel('Glucosio [mg/dL]');
title(['Confronto CGM: ML vs NL - Paziente ' num2str(patient)]);
legend('show');
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white');

%------------IOB------------
%calcolo IOB_cap
IOB_cap_ML = theta_ott_ML(5)*(x_ML(2,:) + x_ML(3,:));
%IOB_cap_NL = theta_ott_NL(4)*(x_NL(2,:) + x_NL(3,:));
IOB_cap_NL = IOB_NL; % già calcolato nella sim. del sistema per deltaIOB

% Plot dei dati 
subplot(3, 1, 2); 
plot(1:1:time, IOB(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'IOB Reale');
hold on;
plot(1:1:time, IOB_cap_NL(1:time), 'k-', 'LineWidth', 1.5, 'DisplayName', 'IOB NL');
plot(1:1:time, IOB_cap_ML(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'IOB ML');

hold off;
grid on;
xlim([0, time]);
xlabel('Tempo [min]');
ylabel('IOB [U]');
title(['Confronto IOB: ML vs NL - Paziente ' num2str(patient)]);
legend('show');
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white');

%-----------Ra-------------
%calcolo Ra_cap
Ra_cap_ML = theta_ott_ML(4)*x_ML(4,:);
Ra_cap_NL = theta_ott_NL(3)*x_NL(4,:);

% Plot dei dati 
subplot(3, 1, 3); 
plot(1:1:time, Ra(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Ra Reale');
hold on;
plot(1:1:time, Ra_cap_NL(1:time), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ra NL');
plot(1:1:time, Ra_cap_ML(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Ra ML');

hold off;
grid on;
xlim([0, time]);
xlabel('Tempo [min]');
ylabel('Ra [mg/(dl\cdotmin)]');
title(['Confronto Ra: ML vs NL - Paziente ' num2str(patient)]);
legend('show');
set(gca, 'FontSize', 12);
set(gcf, 'Color', 'white');

end

