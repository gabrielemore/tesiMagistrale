function [theta_ott_NL] = model_identification_NL(data_circadian_on,patient,div_train_test_gg,theta_ott_ML,grafici,iter)
%% ESTRAZIONE DATI
[time,Gb,~,CF,Ub,u,r,y,IOB,Ra] = data_extraction(data_circadian_on,patient);

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
T_train=60*24*div_train_test_gg;

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
if iter
    options.Display='iter-detailed';
else
    options.Display='none';
end

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
deltaG = zeros(time,1);
deltaIOB = zeros(time,1);
IOB_Sim = zeros(time,1);

i = 1:time;
%baseline insulin sensitivity
Si_tar = x0(6)*(1+0.01*o8*sin((2*pi*i*Ts)/(60*24) + 2*pi*0.01*o9)); % x0(6)=CF
%non abbiamo G_b fisso ma variabile
G_basal = (o0-Si_tar*x0(2))/o1; %x0(2)=Ub;

for k=1:time
    %deviation of gluscose from its basal
    deltaG(k) = x(1,k) - G_basal(k);

    %deviation of IOB from its basal
    IOB_Sim(k) = o3*(x(2,k) + x(3,k));
    deltaIOB(k) = IOB_Sim(k) - IOB_basal;

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
disp("---GoF MODELLO NON LINEARE---");
GoF_train = 100*(1 - (norm(y(1:T_train)- y_cap(1:T_train)))/(norm(y(1:T_train) - mean(y(1:T_train)))));
GoF_test = 100*(1 - (norm(y(T_train+1:time)- y_cap(T_train+1:time)))/(norm(y(T_train+1:time) - mean(y(T_train+1:time)))));
disp("GoF train: " + GoF_train );
disp("GoF test: " + GoF_test );

%% CONFRONTO GRAFICO
if grafici
    %-------------CGM--------------
    figure('Name', ['Confronto CGM: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
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
    title(['Confronto CGM: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');

    %export_fig('identification_NL_CIRC_ON_GA_WORST.pdf','-pdf');

    %------------IOB------------
    %calcolo IOB_cap
    IOB_cap = o3*(x(2,:) + x(3,:));

    figure('Name', ['Confronto IOB: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
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
    title(['Confronto IOB: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');

    %export_fig('identification_NL_CIRC_ON_IOB_WORST.pdf','-pdf');

    %-----------Ra-------------
    %calcolo Ra_cap
    Ra_cap = o2*x(4,:);

    figure('Name', ['Confronto Ra: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
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
    title(['Confronto Ra: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');

    %export_fig('identification_NL_CIRC_ON_RA_WORST.pdf','-pdf');
end

end

