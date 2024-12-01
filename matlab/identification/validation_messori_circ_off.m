function validation_messori_circ_off(data_circadian_off_MESSORI,patient,theta_ott_ML,grafici)
%% ESTRAZIONE DATI
[time,Gb,~,~,Ub,u,r,y,IOB,Ra] = data_extraction(data_circadian_off_MESSORI,patient);

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

%% SIMULAZIONE SISTEMA LINEARE CON PARAMETRI STIMATI
%stati iniziali
x0=[Gb Ub Ub 0 0]';
%definizione tempo di campionamento
Ts = 1;
y_cap = zeros(time,1);
x(:,1)=x0;
for k=1:time
    %sistema lineare
    x(:,k+1) = x(:,k) + Ts*(A*x(:,k) + B_u*u(k) + B_r*r(k) + E);
    y_cap(k)=C*x(:,k);
end

%% METRICHE GoF
disp("---GoF MESSORI VAL. CIRC. OFF---");
GoF = 100*(1 - (norm(y(1:time)- y_cap(1:time)))/(norm(y(1:time) - mean(y(1:time)))));
disp("GoF: " + GoF );

%% CONFRONTO GRAFICO
if grafici
    %-------------CGM--------------
    figure('Name', ['Validazione CGM: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
    % Plot dei dati
    plot(1:1:time, y(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'CGM Reale');
    hold on;
    plot(1:1:time, y_cap(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'CGM Stimato');

    hold off;
    grid on;
    xlim([0, time]);

    xlabel('Tempo [min]');
    ylabel('Glucosio [mg/dL]');
    title(['Validazione CGM: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');

    %------------IOB------------
    %calcolo IOB_cap
    IOB_cap = o4*(x(2,:) + x(3,:));

    figure('Name', ['Validazione IOB: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
    % Plot dei dati
    plot(1:1:time, IOB(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'IOB Reale');
    hold on;
    plot(1:1:time, IOB_cap(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'IOB Stimato');

    hold off;
    grid on;
    xlim([0, time]);
    xlabel('Tempo [min]');
    ylabel('IOB [U]');
    title(['Validazione IOB: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');

    %-----------Ra-------------
    %calcolo Ra_cap
    Ra_cap = o3*x(4,:);

    figure('Name', ['Validazione Ra: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
    % Plot dei dati
    plot(1:1:time, Ra(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Ra Reale');
    hold on;
    plot(1:1:time, Ra_cap(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Ra Stimato');

    hold off;
    grid on;
    xlim([0, time]);
    xlabel('Tempo [min]');
    ylabel('Ra [mg/(dl\cdotmin)]');
    title(['Validazione Ra: Reale vs Stimato - ML - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');
end
end

