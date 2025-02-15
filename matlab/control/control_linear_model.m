function control_linear_model(data_circadian,patient,theta_ott_ML,Ts,deltaT)
%% IMPORT CASADI
import casadi.*

%% ESTRAZIONE DATI
[~,Gb,~,~,Ub,~,~,~,~,~] = data_extraction(data_circadian,patient);

%% DEFINIZIONE SISTEMA LINEARE
o0=theta_ott_ML(1);
o1=theta_ott_ML(2);
o2=theta_ott_ML(3);
o3=theta_ott_ML(4);
o4=theta_ott_ML(5);
o5=theta_ott_ML(6);

%CR fisiologico
CR=o2/o3;

%Forma matriciale
A = [-o1 -o2 0 o3 0;0 -1/o4 1/o4 0 0;0 0 -1/o4 0 0; 0 0 0 -1/o5 1/o5; 0 0 0 0 -1/o5];
B_u= [0 0 1/o4 0 0]';
B_r= [0 0 0 0 1/o5]';
E=[o0 0 0 0 0]';
C=[1 0 0 0 0];

%% DISCRETIZZAZIONE SISTEMA LINEARE
[A_d,B_u_d,B_r_d,C_d,E_d] = system_discretization(A,B_u,B_r,E,C,Ts,deltaT);

%% COSTANTI MPC
%orizzionte predittivo
N=72;

%vincoli controllo
umin = 0;
umax = 15;

%vincoli (nb:x4 e x5 senza vincoli)
x1min = 54;
x1max = 300;
% IOB
CHO_UB = 90;
tau = 120;
IOB_s = o4*2*Ub; %[22 - 6)
IOB_d = IOB_s + (CHO_UB/CR + tau*Ub); %[6 - 22)

%uscita
ymin = 70;
ymax = 140;

%Pesi ipoglicemia e iperglicemia
p_ipo = 1*10^9;
p_iper = 1*10^7;
%Pesi variabili slack stati x1 e IOB(x2 e x3)
p_ipo_x1 =1*10^9;
p_iper_x1 =1*10^9;
p_IOB = 1*10^9;

%% VARIABILI e PARAMETRI OTTIMIZZAZIONE
opti = casadi.Opti();

%variabili di ottimizzazione
xf = opti.variable(5,N+1);
uf = opti.variable(N);

%var. ausiliarie
xa = opti.variable(3,1);
ua = opti.variable();

%costo addizionale variabili slack
delta_iper_x1 = opti.variable(N);
delta_ipo_x1 = opti.variable(N);
delta_IOB = opti.variable(N);

%costo addizionale
delta_ipo = opti.variable();
delta_iper = opti.variable();

%vincolo dinamico IOB(t)
IOB_bound = opti.parameter(N+1);
xk = opti.parameter(5,1); % x0 stimato da osservatore
rk = opti.parameter(N);

%% FUNZIONE DI COSTO
Q=1;
R=100;

V=0;
for j=1:N
    %funzione costo
    V_dyn = (C_d*xf(:,j) - xa(1))'*Q*(C_d*xf(:,j) - xa(1)) + (uf(j) - ua)'*R*(uf(j) - ua);
    %costo slack IOB
    V_IOB = p_IOB*delta_IOB(j)^2;
    %costo slack x1
    V_x1 = p_iper_x1*delta_iper_x1(j)^2 + p_ipo_x1*delta_ipo_x1(j)^2;
    %costo tot
    V = V + V_dyn + V_IOB + V_x1;
end
%costo addizionale
V_s = p_iper*delta_iper^2 + p_ipo*delta_ipo^2;
%costo totale
V = V + V_s;

opti.minimize(V);

%% VINCOLI

%condizione iniziale x0
opti.subject_to(xf(:,1)==xk);

for j=1:N
    %dinamica modello
    opti.subject_to(xf(:,j+1)==A_d*xf(:,j) + B_u_d*uf(j) + B_r_d*rk(j) + E_d);
    %x1 +- slack
    opti.subject_to(xf(1,j) >=x1min - delta_ipo_x1(j));
    opti.subject_to(xf(1,j) <=x1max + delta_iper_x1(j));
    %x2 e x3 +- slack
    %opti.subject_to(o4*(xf(2,j)+xf(3,j)) >=0);
    opti.subject_to(o4*(xf(2,j)+xf(3,j)) <=IOB_bound(j) + delta_IOB(j));
    %u
    opti.subject_to(uf(j) <=umax);
    opti.subject_to(uf(j) >=umin);
    %slack x1 e IOB >=0
    opti.subject_to(delta_ipo_x1(j) >=0);
    opti.subject_to(delta_iper_x1(j) >=0);
    opti.subject_to(delta_IOB(j) >=0);
end

%eq
%x1a = ya
opti.subject_to(xa(1) >=ymin - delta_ipo);
opti.subject_to(xa(1) <=ymax + delta_iper);
%delta ipo e iper >=0
opti.subject_to(delta_ipo >=0);
opti.subject_to(delta_iper >=0);
%x2 e x3
%opti.subject_to(o4*(xa(2)+xa(3)) >=0);
opti.subject_to(o4*(xa(2) + xa(3)) <= IOB_bound(N+1));
%ua
opti.subject_to(ua <=umax);
opti.subject_to(ua >=umin);
%C_n*x(N)=xa
opti.subject_to(xf(1:3,N+1)==xa);
opti.subject_to(xa==A_d(1:3,1:3)*xa + B_u_d(1:3)*ua + E_d(1:3));

%% SOLVER

% Solver options
plgopt = struct;
s_opts = struct;
s_opts.tol = 1e-2;
s_opts.constr_viol_tol = 1e-2;
s_opts.acceptable_tol = 1e-2;
s_opts.acceptable_constr_viol_tol = 0.01;
s_opts.max_iter = 1e4;
s_opts.print_level = 0;
opti.solver('ipopt', plgopt, s_opts);

%funzione MPC
%opti.to_function(name, inputs, outputs, input_names, output_names)
MPC=opti.to_function('MPC',{xk,rk,IOB_bound},{uf,xf,V,delta_ipo,delta_iper,delta_ipo_x1,delta_iper_x1,delta_IOB,xa,ua},{'xk','rk','IOB_bound'},{'uf_opt','xf_opt','V_opt','delta_ipo','delta_iper','delta_ipo_x1','delta_iper_x1','delta_IOB','xa','ua'});

%% OSSERVATORE
%KF settings
ki=1*10^-2;
kd=1;
%varianza
sigma_q = 1;
sigma_q_cgm = 6.51;

% Matrici sistema aumentato
B_0 = [1 0 0 0 0]';
Bd=[B_u_d B_r_d E_d];

Aaug = [A_d B_0 zeros(size(A_d,1),1);
    zeros(1,size(A_d,1)), 1, Ts;
    zeros(1,size(A_d,1)), 0, 1];
Baug = [Bd;
    zeros(2,size(Bd,2))];
G = [0 ki ki 0 0 0 kd]';
Caug = [C_d 0 0];

Q_kf = G * sigma_q * G';
R_kf = sigma_q_cgm;
%P(0|0)=P(inf|inf)
[~,Pkk,~,~] = dlqe(Aaug,G,Caug,sigma_q,R_kf);

%% SIMULZIONE
%vettori per salvataggio dati intermedi
v_xf =[];
v_xa =[];
v_uf =[];
v_ua =[];
v_u=[];
v_y=[];
v_exit=[];
v_VN=[];
v_x_sim=[];
v_tt =[];
v_pkk = [];
v_kk = [];
v_d1 = [];
v_d2 =[];
v_delta_ipo = [];
v_delta_iper = [];
v_delta_iper_x1 =  [];
v_delta_ipo_x1 =  [];
v_delta_IOB =[];
v_xk_obs =[];
v_x_real =[];
v_time_sol =[];

%tempo di simulazione in minuti
Tmax = 24*60*1;

%ipotesi iniziali
x0 = [Gb Ub Ub 0 0]';
u0 = 0;

%costruzione upper-lower bound IOB
IOB_vet = create_IOB_vector(Tmax,IOB_s,IOB_d);

%simulazione pasti (1 giorno + 6h)
% rk_in = zeros(1,Tmax+(N*Ts));
[rk_in,~] = create_RK_random(Tmax+(N*Ts));
%rk_in(801:830) = 60/30; %60g divisi in 30 minuti
%rk_in(1001:1020) = 30/20; %30g divisi in 20 minuti

%insulina da iniettare
u_in = 0;
%intervallo simulazione ode45
tspan = [0 1];

%inizializzazioni osservatore
x_aug = [x0;0;0];
u_aug = [u0;0;1];
y=Gb;
%inizializzazione primo stato sistema reale
x_real = x0;

for k=1:Tmax
    if mod(k-1,Ts) == 0
        %OSSERVATORE
        %stima x0 con osservatore
        [xk_obs,Pkk,Kk]=ODO(Aaug,Baug,Caug,y,x_aug,u_aug,Pkk,Q_kf,R_kf);
        %non considero i disturbi per MPC (prendo solo i primi 5 stati)
        xk_sim = xk_obs(1:5);

        %costruzione vettore rk_sim di lunghezza N (predizione 6h)
        rk_sim = create_RK_Sim(rk_in,k,Ts,N);

        %MPC
        %risolvo il problema di minimizzazione al passo k
        [uf_sol,xf_sol,V_sol,d_ipo,d_iper,d_ipo_x1,d_iper_x1,d_IOB,xa_sol,ua_sol]=MPC(xk_sim,rk_sim,IOB_vet(k:Ts:k+N*Ts));
        %verifico exit MPC
        exit=get_stats(MPC);

        %converto risultato simbolico in risultato numerico
        xf_sim=full(xf_sol);
        uf_sim=full(uf_sol);
        xa_sim=full(xa_sol);
        ua_sim=full(ua_sol);
        Vn=full(V_sol);
        delta_iper_sim = full(d_iper);
        delta_ipo_sim = full(d_ipo);
        delta_iper_x1_sim = full(d_iper_x1);
        delta_ipo_x1_sim = full(d_ipo_x1);
        delta_IOB_sim = full(d_IOB);

        if exit.success==1
            % prendo solo il primo elemento di uf ottima e lo applico al sistema
            u_in=uf_sim(1);
        else
            %se non trovo una soluzione ottima utilizzo la soluzione al
            %passo precedente
            u_in=u_aug(1);
        end


        %matrice stati x controllati con MPC su N=72 (6 ore)
        v_xf = [v_xf;xf_sim];
        %matrice controllo u
        v_uf = [v_uf;uf_sim];
        %matrice stati artificiali
        v_xa = [v_xa;xa_sim];
        %matrice controllo artificiale
        v_ua = [v_ua;ua_sim];
        %matrice valore del costo
        v_VN = [v_VN; Vn];
        %matrice soluzione MPC trovata
        v_exit=[v_exit,exit.success];
        %salvo matrici osservatore per vedere se convergono
        v_kk = [v_kk;Kk];
        v_pkk = [v_pkk;Pkk];
        %salvo valori disturbi per vedere se convergono
        v_d1 = [v_d1;xk_obs(6)];
        v_d2 = [v_d2;xk_obs(7)];
        %salvo valori delta calcolati in funzione MPC
        v_delta_iper = [v_delta_iper;delta_iper_sim];
        v_delta_ipo = [v_delta_ipo;delta_ipo_sim];
        v_delta_iper_x1 = [v_delta_iper_x1;delta_iper_x1_sim];
        v_delta_ipo_x1 = [v_delta_ipo_x1;delta_ipo_x1_sim];
        v_delta_IOB = [v_delta_IOB;delta_IOB_sim];
        %salvo valori xk osservati dall'osservatore
        v_xk_obs = [v_xk_obs xk_obs];
        %salvo tempi risoluzione MPC di ogni iterazione
        v_time_sol = [v_time_sol;exit.t_wall_total];

        %aggiorno valori aug per prossima iterazione MPC
        x_aug = xk_obs;
        u_aug = [u_in;rk_in(k);1];
    end

    %avanzamento sistema reale
    [tt,x_sim] = ode45(@(t, x) patient_fun(t, x, A, B_u, B_r, E, u_in, rk_in(k)), tspan, x_real);
    % Prendo l'ultimo stato
    x_real = x_sim(end, :)';
    %unica misura che sono in grado di effettuare è x1
    y=x_real(1);

    %salvo i valori insulina inserita
    v_u= [v_u,u_in];

    %azzero ingresso insulina per ciclo sucessivo
    u_in = 0;
    %salvo in un vettore le x ottenute dalla simulazione del sistema reale
    v_x_sim=[v_x_sim;x_sim(1:end-1,:)];
    %salvo in un vettore gli istanti temporali corrispondenti alle x della
    %sim
    v_tt=[v_tt; tt(1:end-1) + (k-1)];
    %vettore uscite y
    v_y = [v_y;y];
    %vettore x real
    v_x_real = [v_x_real x_real];
end

%% TEST - salvattaggio vettori per grafici
% Salva i vettori in un file MAT
save('dati_intermedi.mat', 'v_xf', 'v_u', 'v_y', 'v_exit', 'v_VN', ...
    'v_x_sim', 'v_tt', 'v_pkk', 'v_kk', 'v_d1', 'v_d2', ...
    'v_delta_ipo', 'v_delta_iper', 'v_xk_obs','v_x_real','IOB_vet','rk_in');

%% GRAFICI
figure('Name', ['Controllo MPC - Paziente ' num2str(patient)]);

%---------GLICEMIA---------
subplot(5, 1, 1);
plot(v_y(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Glicemia reale');
hold on;
plot(v_xk_obs(1,:), 'b-', 'LineWidth', 1, 'DisplayName', 'Glicemia osservatore (ODO)');
hold off;
grid on;
xlim([1, length(v_xk_obs)]);
ylabel('Glicemia [mg/dL]');
legend('show');
set(gca, 'FontSize', 12);

%------------IOB------------
IOB_obs = o4*(v_xk_obs(2,:) + v_xk_obs(3,:));
IOB_real = o4*(v_x_real(2,:) + v_x_real(3,:));
subplot(5, 1, 2);
yyaxis left;
plot(IOB_vet(1:Ts:end), 'k--', 'LineWidth', 1, 'DisplayName', 'Vincolo IOB');
hold on;
plot(IOB_real(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'IOB reale');
plot(IOB_obs, 'b-', 'LineWidth', 1, 'DisplayName', 'IOB osservatore (ODO)');
ylabel('IOB [U]');
yyaxis right;
stem(v_u(1:Ts:end), 'g', 'filled', 'LineWidth', 1, 'DisplayName', 'Insulina iniettata', 'Marker', 'none');
ylabel('Insulina [U/min]');
hold off;
grid on;
xlim([1, length(IOB_obs)]);
legend('show');
set(gca, 'FontSize', 12);

%-----------Ra-------------
Ra_real = o3*v_x_real(4,:);
Ra_obs = o3*v_xk_obs(4,:);
subplot(5, 1, 3);
yyaxis left;
plot(Ra_real(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Ra reale');
hold on;
plot(Ra_obs, 'b-', 'LineWidth', 1, 'DisplayName', 'Ra osservatore (ODO)');
ylabel('Ra [mg/(dl\cdotmin)]');
yyaxis right;
plot(rk_in(1:Ts:end), 'k-', 'LineWidth', 1, 'DisplayName', 'Pasti');
ylabel('Ingest [g/min]');
hold off;
grid on;
xlim([0, length(v_xk_obs)]);
ylim([0, 5]);
legend('show');
set(gca, 'FontSize', 12);

%-----------D1-------------
subplot(5, 1, 4);
plot(v_d1, 'b-', 'LineWidth', 1, 'DisplayName', 'Disturbo');
grid on;
xlim([0, length(v_d1)]);
ylabel('d_1');
legend('show');
set(gca, 'FontSize', 12);

%---------COSTO----------
subplot(5, 1, 5);
plot(v_VN, 'b-', 'LineWidth', 1, 'DisplayName', 'Costo');
grid on;
xlim([0, length(v_VN)]);
%Etichetta asse X in comune
xlabel('Istanti [min]');
ylabel('Vn');
legend('show');
set(gca, 'FontSize', 12);

%Sfondo bianco per l'intera figura
set(gcf, 'Color', 'white');
