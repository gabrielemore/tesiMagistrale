function control_linear_model(data_circadian,patient,theta_ott_ML,Ts,deltaT)
%% IMPORT CASADI
import casadi.*

%% ESTRAZIONE DATI
[~,Gb,CR,~,Ub,~,~,~,~,~] = data_extraction(data_circadian,patient);

%% DEFINIZIONE SISTEMA LINEARE
o0=theta_ott_ML(1);
o1=theta_ott_ML(2);
o2=theta_ott_ML(3);
o3=theta_ott_ML(4);
o4=theta_ott_ML(5);
o5=theta_ott_ML(6);

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
p_ipo = 1*10^7;
p_iper = 1*10^9;

%% VARIABILI e PARAMETRI OTTIMIZZAZIONE
opti = casadi.Opti();

%variabili di ottimizzazione
xf = opti.variable(5,N+1);
uf = opti.variable(1,N);

%var. ausiliarie
xa = opti.variable(3,1);
ua = opti.variable(1,1);

%costo addizionale
delta_ipo = opti.variable(1,1);
delta_iper = opti.variable(1,1);

%vincolo dinamico IOB(t)
IOB_bound = opti.parameter(1,N+1);
xk = opti.parameter(5,1); % x0 stimato da osservatore
rk = opti.parameter(1,N);

%% FUNZIONE DI COSTO
Q=eye(1);
R=eye(1);

V=0;
for j=1:N
    V_dyn = (C_d*xf(:,j) - xa(1))'*Q*(C_d*xf(:,j) - xa(1)) + (uf(j) - ua)'*R*(uf(j) - ua);
    V = V + V_dyn;
end
%costo addizionale
V_s = p_iper*delta_iper^2 + p_ipo*delta_ipo^2;
V = V + V_s;

opti.minimize(V);

%% VINCOLI

%condizione iniziale x0
opti.subject_to(xf(:,1)==xk);

for j=1:N
    %dinamica modello
    opti.subject_to(xf(:,j+1)==A_d*xf(:,j) + B_u_d*uf(j) + B_r_d*rk(j) + E_d);
    %x1
    opti.subject_to(xf(1,j) >=x1min);
    opti.subject_to(xf(1,j) <=x1max);
    %x2 e x3
    opti.subject_to(o4*(xf(2,j)+xf(3,j)) >=0);
    opti.subject_to(o4*(xf(2,j) + xf(3,j)) <= IOB_bound(j));
    %u
    opti.subject_to(uf(j) <=umax);
    opti.subject_to(uf(j) >=umin);
end

%eq
%x1a = ya
opti.subject_to(xa(1) >=ymin - delta_ipo);
opti.subject_to(xa(1) <=ymax + delta_iper);
%delta ipo e iper >=0
opti.subject_to(delta_ipo >= 0);
opti.subject_to(delta_iper >=0);
%x2 e x3
opti.subject_to(o4*(xa(2)+xa(3)) >=0);
opti.subject_to(o4*(xa(2) + xa(3)) <= IOB_bound(N+1));
%ua
opti.subject_to(ua <=umax);
opti.subject_to(ua >=umin);
%C_n*x(N)=xa
opti.subject_to(xf(1:3,N+1)==xa);
opti.subject_to(xa==A_d(1:3,1:3)*xa + B_u_d(1:3)*ua + E_d(1:3));

%% SOLVER
opti.solver('ipopt');

%funzione MPC
%opti.to_function(name, inputs, outputs, input_names, output_names)
MPC=opti.to_function('MPC',{xk,rk,IOB_bound},{uf,xf,V},{'xk','rk','IOB_bound'},{'uf_opt','xf_opt','V_opt'});

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
%[~,Pkk,~,~] = dlqe(Aaug,G,Caug,Q_kf,R_kf); ???
[~,Pkk,~,~] = dlqe(Aaug,G,Caug,sigma_q,R_kf);

%% SIMULZIONE
% vx=[];
% vu=[];
% vys=[];
% vexit=[];
% vVN=[];
% Vn=0;
vx_sim=[];

%tempo di simulazione (1 giorno)
Tmax = 1440;

%condizioni iniziali
x0 = [Gb Ub Ub 0 0]';
u0 = 0;

%costruzione upper-lower bound IOB + 6h
IOB_vet = zeros(1,1800);
IOB_vet(1:360) = IOB_s;
IOB_vet(361:1320) = IOB_d;
IOB_vet(1321:1440) = IOB_s;
% + 6h
IOB_vet(1441:1800) = IOB_s;

%simulazione pasti
rk = zeros(1,Tmax);%da corregere, dovrebbe essere fino a 1440 e poi mando di volta in volta 72 elementi a MPC

%osservatore
x_aug = [x0;0;0];
u_aug = [u0;rk(1);1];
y=x0(1);

%insulina da iniettare
u_in = 0;

%intervallo simulazione ode45
tspan = [0 1];
%inizializzo primo stato sistema
x_real = x0;
 
% vys=[vys,xk(1)];
% vx=[vx,xk];

for k=1:Tmax

    if mod(k,Ts) == 0
        %stima x0 con osservatore
        [xk_obs,Pkk]=ODO(Aaug,Baug,Caug,y,x_aug,u_aug,Pkk,Q_kf,R_kf);
        %non considero i disturbi per MPC
        xk = xk_obs(1:5);
        %risolvo il problema di minimizzazione al passo k
        [uf_sol,xf_sol,V_sol]=MPC(xk,rk(k:Ts:k+N*Ts-1),IOB_vet(k:Ts:k+N*Ts));

        % exit=get_stats(MPC);
        % vexit=[vexit,exit.success];

        %converto risultato simbolico in risultato numerico
        uf=full(uf_sol);
        % xf=full(xf_sol);
        % Vn=full(V_sol);

        % prendo solo il primo elemento di uf ottima e lo applico al sistema
        u_in=uf(:,1);
        %valori aug per prossima iterazione aggiornati
        x_aug = xk_obs;
        u_aug = [u_in;rk(k);1];
    end
    
    %avanzamento sistema reale
    [tt,x_sim] = ode45(@(t, x) patient_fun(t, x, A, B_u, B_r, E, u_in, rk(k)), tspan, x_real);
    % Prendo l'ultimo stato
    x_real = x_sim(end, :)'; 
    %unica misura che sono in grado di effettuare è x1
    y=x_real(1);
    %azzero ingresso per ciclo sucessivo
    u_in = 0;
    %salvo in un vettore x ottenuti dalla simulazione del sistema reale
    vx_sim=[vx_sim;x_sim];
end

end

