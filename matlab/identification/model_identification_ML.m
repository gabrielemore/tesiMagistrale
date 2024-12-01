function [theta_ott_ML] = model_identification_ML(data_circadian_off,patient,o4_ini,div_train_test_gg,grafici,iter)
%% ESTRAZIONE DATI
[time,Gb,CR,CF,Ub,u,r,y,IOB,Ra] = data_extraction(data_circadian_off,patient);

%% STIMA DEI PARAMETRI
%inizializzazione valori iniziali parametri
o01 = 0.004; % [1/min]
o02 = CF; % [mg/(dL·U)]
o03 = CF/CR; % [mg/(dL·g)]
o04 = o4_ini; % [min]
o05 = 40; % [min]
o00 = o01*Gb + o02*Ub; % [mg/dL/min]

%vettore teta0 iniziale
theta0= [o00 o01 o02 o03 o04 o05]';
%divisione dataset train e test
T_train=60*24*div_train_test_gg;

%definizione tempo di campionamento
Ts = 1;
%stati iniziali
x0=[Gb Ub Ub 0 0]';

min = zeros(6,1);
max = ones(6,1)*inf;

%options = optimoptions('fmincon','Algorithm','interior-point','EnableFeasibilityMode',true, SubproblemAlgorithm='cg');
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxIterations=1000;
options.FunctionTolerance = 1e-5;
options.StepTolerance = 1e-10;
options.MaxFunctionEvaluations=100000;
if iter
    options.Display='iter-detailed';
else
    options.Display='none';
end
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
disp("---PARAMETRI MODELLO LINEARE---");
disp("CR: " + CR + " - CR stimato: " + CR_cap);
disp("theta4 iniziale: " + o04 + " - theta4 stimato: " + o4);
disp("theta5 iniziale: " + o05 + " - theta5 stimato: " + o5);

%% SIMULAZIONE SISTEMA LINEARE CON PARAMETRI STIMATI
y_cap = zeros(time,1);
x(:,1)=x0;
for k=1:time
    %sistema lineare
    x(:,k+1) = x(:,k) + Ts*(A*x(:,k) + B_u*u(k) + B_r*r(k) + E);
    y_cap(k)=C*x(:,k);
end

%% METRICHE GoF TRAIN E TEST
disp("---GoF MODELLO LINEARE---");
GoF_train = 100*(1 - (norm(y(1:T_train)- y_cap(1:T_train)))/(norm(y(1:T_train) - mean(y(1:T_train)))));
GoF_test = 100*(1 - (norm(y(T_train+1:time)- y_cap(T_train+1:time)))/(norm(y(T_train+1:time) - mean(y(T_train+1:time)))));
disp("GoF train: " + GoF_train );
disp("GoF test: " + GoF_test );

%% CONFRONTO GRAFICO
if grafici
    %-------------CGM--------------
    figure('Name', ['Confronto CGM: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
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
    title(['Confronto CGM: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');

    %------------IOB------------
    %calcolo IOB_cap
    IOB_cap = o4*(x(2,:) + x(3,:));

    figure('Name', ['Confronto IOB: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
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
    title(['Confronto IOB: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');

    %-----------Ra-------------
    %calcolo Ra_cap
    Ra_cap = o3*x(4,:);

    figure('Name', ['Confronto Ra: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
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
    title(['Confronto Ra: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');
end

end

