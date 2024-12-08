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

C_n = [eye(3) zeros(3,2)];

%% COSTANTI MPC
%orizzionte predittivo
N=72;

%vincoli controllo
umin = 0;
umax = 15;

%vincoli stati (x4 e x5 senza vincoli)
x1min = 60;
x1max = 300;

% ????????????? da verificare
CHO_UB = 80;
tau = 5;
% ????????????? da verificare

IOB_s = o4*2*Ub; %[6 - 22)
IOB_d = IOB_s + (CHO_UB/CR + tau*Ub); %[22 - 6)

%uscita
ymin = 70;
ymax = 140;

%?????? Pesi ipoglicemia e iperglicemia
p_ipo = 1;
p_iper = 0.1;

%% VARIABILI SIMBOLICHE CasADi
x1 = MX.sym('x1');
x2 = MX.sym('x2');
x3 = MX.sym('x3');
x4 = MX.sym('x4');
x5 = MX.sym('x5');

x=[x1;x2;x3;x4;x5];

u = MX.sym('u');
r = MX.sym('r');

%% DINAMICA SISTEMA
x_next = A_d*x + B_u_d*u + B_r_d*r + E_d;
% (nome simb. fun. - var. input - var. output - nomi sim. var - nome simb.
% out)
sis=Function('sis',{x,u,r},{x_next},{'x','u','r'},{'x_next'});

%% VARIABILI e PARAMETRI OTTIMIZZAZIONE
opti = casadi.Opti();

%variabili di ottimizzazione (vedere se modificare 5 e 1 con variabili
%prese da matrice B)
xf = opti.variable(5,N+1);
uf = opti.variable(1,N);

%equilibrio
xa = opti.variable(5,1);
ua = opti.variable(1,1);

%costo addizionale
delta_ipo = opti.variable(1,1);
delta_iper = opti.variable(1,1);

%t_attuale necessario per vincolo dinamico IOB(t)
t_attuale = opti.parameter(1,1);
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

V_s = p_iper*delta_iper^2 + p_ipo*delta_ipo^2;
V = V + V_s;

opti.minimize(V);

%% VINCOLI

%condizione iniziale x0
opti.subject_to(xf(:,1)==xk);

for j=1:N
    %dinamica modello
    opti.subject_to(xf(:,j+1)==sis(xf(:,j),uf(j),rk(j)));
    %x1
    opti.subject_to(xf(1,j) >=x1min);
    opti.subject_to(xf(1,j) <=x1max);
    %x2 e x3
    opti.subject_to(o4*(xf(2,j)+xf(3,j)) >=0);
    value = fmod(j + t_attuale, 1440);
    %opti.subject_to(o4*(xf(2,j) + xf(3,j)) <= IOB_d * (360 <= value <= 1320) + IOB_s * ~(360 <= value <= 1320));
    opti.subject_to(o4*(xf(2,j) + xf(3,j)) <= IOB_d * (value >= 360 & value <= 1320) + IOB_s * (value < 360 | value > 1320));
    %x4 e x5
    opti.subject_to(xf(4:5,j) >=0);
    %u
    opti.subject_to(uf(j) <=umax);
    opti.subject_to(uf(j) >=umin);
end

%var artificiale equilibrio
%x1a = ya
opti.subject_to(xa(1) >=ymin - delta_ipo);
opti.subject_to(xa(1) <=ymax + delta_iper);
%delta ipo e iper >=0
opti.subject_to(delta_ipo >= 0);
opti.subject_to(delta_iper >=0);
%x2 e x3
opti.subject_to(o4*(xa(2)+xa(3)) >=0);
value = fmod(j + t_attuale, 1440);
%opti.subject_to(o4*(xa(2) + xa(3)) <= IOB_d * (360 <= value <= 1320) + IOB_s * ~(360 <= value <= 1320));
opti.subject_to(o4*(xa(2) + xa(3)) <= IOB_d * (value >= 360 & value <= 1320) + IOB_s * (value < 360 | value > 1320));
%x4 e x5
opti.subject_to(xa(4:5) >=0);
%ua
opti.subject_to(ua <=umax);
opti.subject_to(ua >=umin);
%C_n*x(N)=xa
opti.subject_to(C_n*xf(:,N+1)==xa);
opti.subject_to(xa==sis(xa,ua,0));

%% SOLVER
opti.solver('ipopt');

%funzione MPC
%opti.to_function(name, inputs, outputs, input_names, output_names)
MPC=opti.to_function('MPC',{xk,rk,t_attuale},{uf,xf,V},{'xk','rk','t_attuale'},{'uf_opt','xf_opt','V_opt'});

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

Q_kf = G * sigma_q * G';
R_kf = sigma_q_cgm;
%P(0|0)=P(inf|inf)
[Pkk, ~, ~] = idare(Aaug, Baug, Q_kf, R_kf,[],[]);

%% SIMULZIONE
vx=[];
vu=[];
vys=[];
vexit=[];
vVN=[];
Vn=0;


Tmax = 1440; % tempo massimo di simulazione (1 giorno)

x0=[Gb Ub Ub 0 0]';

xk=x0;
uk=0;
rk = zeros(1,N);

%osservatore
x_aug = [x0;0;0];

vys=[vys,xk(1)];
vx=[vx,xk];

for k=1:Tmax

    [uf_sol,xf_sol,V_sol]=MPC(xk,rk,k);

    exit=get_stats(MPC);
    vexit=[vexit,exit.success];

    %converto risultato simbolico in risultato numerico
    uf=full(uf_sol);
    xf=full(xf_sol);
    Vn=full(V_sol);

    % prendo solo il primo elemento di uf ottima e lo applico al sistema
    uk=uf(:,1);
    x_sis=full(sis(xk,uk,rk)); % faccio avanzare il sistema

    % stimo lo stato iniziale sucessivo con osservatore
    u_aug = [uk;rk(k);1];
    y=x_sis(1);%unica misura che sono in grado di effettuare è x1
    [x_aug,Pkk] = ODO(Aaug,Baug,C_d,y,x_aug,u_aug,Pkk);
    
    xk=x_aug(1:5);%stato stimato dall'osservatore
    
    vx=[vx,xk];
    vu=[vu,uk];
    vVN=[vVN,Vn];
end


figure
plot(vys(1,:))
hold on
plot(vx(1,:))
figure
plot(vVN)
figure
plot(vu(1,:))
hold on
plot(vu(2,:))

end

