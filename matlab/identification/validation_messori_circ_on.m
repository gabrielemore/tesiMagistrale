function validation_messori_circ_on(data_circadian_on_MESSORI,patient,theta_ott_NL,grafici)
%% ESTRAZIONE DATI
[time,Gb,~,CF,Ub,u,r,y,IOB,Ra] = data_extraction(data_circadian_on_MESSORI,patient);

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
deltaG = zeros(time,1);
deltaIOB = zeros(time,1);
IOB_NL = zeros(time,1);

i = 1:time;
%baseline insulin sensitivity
Si_tar = x0(6)*(1+0.01*o8*sin((2*pi*i*Ts)/(60*24) + 2*pi*0.01*o9)); % x0(6)=CF
%non abbiamo G_b fisso ma variabile
G_basal = (o0-Si_tar*x0(2))/o1; %x0(2)=Ub;

for k=1:time
    %deviation of gluscose from its basal
    deltaG(k) = x(1,k) - G_basal(k);

    %deviation of IOB from its basal
    IOB_NL(k) = o3*(x(2,k) + x(3,k));
    deltaIOB(k) = IOB_NL(k) - IOB_basal;

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
disp("---GoF MESSORI VAL. CIRC. ON---");
GoF = 100*(1 - (norm(y(1:time)- y_cap(1:time)))/(norm(y(1:time) - mean(y(1:time)))));
disp("GoF: " + GoF );

%% CONFRONTO GRAFICO
if grafici
    %-------------CGM--------------
    figure('Name', ['Validazione CGM: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
    % Plot dei dati
    plot(1:1:time, y(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'CGM Reale');
    hold on;
    plot(1:1:time, y_cap(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'CGM Stimato');

    hold off;
    grid on;
    xlim([0, time]);

    xlabel('Tempo [min]');
    ylabel('Glucosio [mg/dL]');
    title(['Validazione CGM: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');

    %------------IOB------------
    %calcolo IOB_cap
    IOB_cap = o3*(x(2,:) + x(3,:));

    figure('Name', ['Validazione IOB: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
    % Plot dei dati
    plot(1:1:time, IOB(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'IOB Reale');
    hold on;
    plot(1:1:time, IOB_cap(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'IOB Stimato');

    hold off;
    grid on;
    xlim([0, time]);
    xlabel('Tempo [min]');
    ylabel('IOB [U]');
    title(['Validazione IOB: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');

    %-----------Ra-------------
    %calcolo Ra_cap
    Ra_cap = o2*x(4,:);

    figure('Name', ['Validazione Ra: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
    % Plot dei dati
    plot(1:1:time, Ra(1:time), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Ra Reale');
    hold on;
    plot(1:1:time, Ra_cap(1:time), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Ra Stimato');

    hold off;
    grid on;
    xlim([0, time]);
    xlabel('Tempo [min]');
    ylabel('Ra [mg/(dl\cdotmin)]');
    title(['Validazione Ra: Reale vs Stimato - NL - Paziente ' num2str(patient)]);
    legend('show');
    set(gca, 'FontSize', 12);
    set(gcf, 'Color', 'white');
end
end

