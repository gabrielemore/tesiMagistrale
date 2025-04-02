function [MPC,N,IOB_s,IOB_d] = initializeProblemMPC_NL(o0,o1,o2,o3,o4,o5,o6,o7,Ub,CR,T,Ts,CF)

addpath('C:\Users\cal_p\Documents\MATLAB\casadi-3.6.5');
import casadi.*

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
%non è piu o4 ma o3
IOB_s = o3*2*Ub; %[22 - 6)
IOB_d = IOB_s + (CHO_UB/CR + tau*Ub); %[6 - 22)

%uscita
% ymin = 90;
ymin = 70;
ymax = 180;

%Pesi ipoglicemia e iperglicemia
p_ipo = 1*10^8;
p_iper = 1*10^6;
%Pesi variabili slack stati x1 e IOB(x2 e x3)
p_ipo_x1 =1*10^8;
p_iper_x1 =1*10^8;
p_IOB = 1*10^8;

p_XN = eye(4)*1*10^4;
p_SIM= eye(4)*1*10^4;





opti = casadi.Opti();

%variabili di ottimizzazione
xf = opti.variable(6,N+1);
uf = opti.variable(N);

%var. ausiliarie
xa = opti.variable(4,T+1);
ua = opti.variable(1,T);

%costo addizionale variabili slack
delta_iper_x1 = opti.variable(N+1);
delta_ipo_x1 = opti.variable(N+1);
delta_IOB = opti.variable(N+1);

%costo addizionale
delta_ipo = opti.variable(T+1);
delta_iper = opti.variable(T+1);
delta_IOB_eq = opti.variable(T+1);

delta_XN = opti.variable(4,1);
delta_simm = opti.variable(4,1);

%vincolo dinamico IOB(t)
IOB_bound = opti.parameter(N+1);
xk = opti.parameter(6,1); % x0 stimato da osservatore
rk = opti.parameter(N);

%parte NON LINEARE
%parametri
Si_Tar_p = opti.parameter(N);
G_basal_p = opti.parameter(N);
IOB_basal_p = opti.parameter();
%parametri per traiettoria equilibrio
Si_Tar_eq = opti.parameter(T);
G_basal_eq = opti.parameter(T+1);
IOB_bound_eq = opti.parameter(T+1);




Q=10^5;
R=10^2;
S=10^5;

V=0;
for j=1:N
    %funzione costo
    V_dyn = (xf(1,j) - xa(1,j))'*Q*(xf(1,j) - xa(1,j)) + (uf(j) - ua(j))'*R*(uf(j) - ua(j));
    %costo slack IOB
    V_IOB = p_IOB*delta_IOB(j)^2;
    %costo slack x1
    V_x1 = p_iper_x1*delta_iper_x1(j)^2 + p_ipo_x1*delta_ipo_x1(j)^2;
    %costo tot
    V = V + V_dyn + V_IOB + V_x1;
end

%costo elementi N+1
V= V + (xf(1,N+1) - xa(1,N+1))'*Q*(xf(1,N+1)- xa(1,N+1)) + p_IOB*delta_IOB(N+1)^2 + p_iper_x1*delta_iper_x1(N+1)^2 + p_ipo_x1*delta_ipo_x1(N+1)^2;


%costo traiettoria di equilibrio
for j=1:T
    V_traj = (xa(1,j) - G_basal_eq(j))'*S*(xa(1,j) - G_basal_eq(j));
    %costo addizionale
    V_s = p_iper*delta_iper(j)^2 + p_ipo*delta_ipo(j)^2;
    %costo slack IOB
    V_IOB_eq = p_IOB*delta_IOB_eq(j)^2;
    %costo totale
    V = V + V_s + V_IOB_eq + V_traj;
end

%costo elemento T+1
V=V + (xa(1,T+1) - G_basal_eq(T+1))'*S*(xa(1,T+1) - G_basal_eq(T+1)) + p_iper*delta_iper(T+1)^2 + p_ipo*delta_ipo(T+1)^2 + p_IOB*delta_IOB_eq(T+1)^2;

V_XN = delta_XN(1:4)'*p_XN*delta_XN(1:4);

V_SIM = delta_simm(1:4)'*p_SIM*delta_simm(1:4);

V = V + V_XN + V_SIM;

opti.minimize(V);

%warm start ua
opti.set_initial(ua, Ub * ones(1, T));
opti.set_initial(xa(1,:), 120 * ones(1, T+1));
opti.set_initial(xa(2,:), Ub * ones(1, T+1));
opti.set_initial(xa(3,:), Ub * ones(1, T+1));
opti.set_initial(xa(4,:), CF * ones(1, T+1));




%condizione iniziale x0
opti.subject_to(xf(:,1)==xk);

for j=1:N
    %dinamica modello non lineare
    opti.subject_to(xf(1,j+1) == xf(1,j) + Ts*(o0 - (o1 * xf(1,j)) - (xf(6,j) * xf(2,j)) + (o2 * xf(4,j)))); % G
    opti.subject_to(xf(2,j+1) == xf(2,j) + Ts*(-(1/o3 * xf(2,j)) + (1/o3 * xf(3,j)))); % Qi
    opti.subject_to(xf(3,j+1) == xf(3,j) + Ts*(-(1/o3 * xf(3,j)) + (1/o3 * uf(j)))); % Qisub
    opti.subject_to(xf(4,j+1) == xf(4,j) + Ts*(-(1/o4 * xf(4,j)) + (1/o4 * xf(5,j)))); % Qg
    opti.subject_to(xf(5,j+1) == xf(5,j) + Ts*(-(1/o4 * xf(5,j)) + (1/o4 * rk(j)))); % Qsto
    opti.subject_to(xf(6,j+1) == xf(6,j) + Ts*(-(1/o5 * xf(6,j)) - (o1/o6 * (xf(1,j) - G_basal_p(j))) - (1/o7 * (o3*(xf(2,j) + xf(3,j)) - IOB_basal_p)) + (1/o5 * Si_Tar_p(j)))); % Si

    %x1 +- slack
    opti.subject_to(xf(1,j) >=x1min - delta_ipo_x1(j));
    opti.subject_to(xf(1,j) <=x1max + delta_iper_x1(j));
    %x2 e x3 +- slack
    opti.subject_to(o3*(xf(2,j)+xf(3,j)) <=IOB_bound(j) + delta_IOB(j));
    %u
    opti.subject_to(uf(j) <=umax);
    opti.subject_to(uf(j) >=umin);
    %slack x1 e IOB >=0
    opti.subject_to(delta_ipo_x1(j) >=0);
    opti.subject_to(delta_iper_x1(j) >=0);
    opti.subject_to(delta_IOB(j) >=0);

    %opti.subject_to(xf(:,j) >=0);
end

%N+1
%x1 +- slack
opti.subject_to(xf(1,N+1) >=x1min - delta_ipo_x1(N+1));
opti.subject_to(xf(1,N+1) <=x1max + delta_iper_x1(N+1));
%x2 e x3 +- slack
opti.subject_to(o3*(xf(2,N+1)+xf(3,N+1)) <=IOB_bound(N+1) + delta_IOB(N+1));
%slack x1 e IOB >=0
opti.subject_to(delta_ipo_x1(N+1) >=0);
opti.subject_to(delta_iper_x1(N+1) >=0);
opti.subject_to(delta_IOB(N+1) >=0);
%opti.subject_to(xf(:,N+1) >=0);


%intersezione traiettoria sistema e traiettoria equilibrio
opti.subject_to((xf(1:3,N+1) - xa(1:3,N+1)) <= delta_XN(1:3));
opti.subject_to(-(xf(1:3,N+1) - xa(1:3,N+1)) <= delta_XN(1:3));

opti.subject_to((xf(6,N+1) - xa(4,N+1)) <= delta_XN(4));
opti.subject_to( -(xf(6,N+1) - xa(4,N+1)) <= delta_XN(4));

opti.subject_to(delta_XN(:)>=0);

%stato iniziale e terminale traiettoria eq devono coincidere
opti.subject_to((xa(:,1) - xa(:,T+1)) <=  delta_simm(:));
opti.subject_to( -(xa(:,1) - xa(:,T+1)) <=  delta_simm(:));

opti.subject_to(delta_simm(:)>=0);

%traiettoria di eq
for j=1:T
    %dinamica modello non lineare all'equilibrio
    opti.subject_to(xa(1,j+1) == xa(1,j) + Ts*(o0 - (o1 * xa(1,j)) - (xa(4,j) * xa(2,j)))); % G
    opti.subject_to(xa(2,j+1) == xa(2,j) + Ts*(-(1/o3 * xa(2,j)) + (1/o3 * xa(3,j)))); % Qi
    opti.subject_to(xa(3,j+1) == xa(3,j) + Ts*(-(1/o3 * xa(3,j)) + (1/o3 * ua(j)))); % Qisub
    opti.subject_to(xa(4,j+1) == xa(4,j) + Ts*(-(1/o5 * xa(4,j)) - (o1/o6 * (xa(1,j) - (o0 - xa(4,j)*ua(j))/o1)) - (1/o7 * (o3*(xa(2,j) + xa(3,j)) - 2*o3*ua(j))) + (1/o5 * Si_Tar_eq(j)))); % Si

    %x1a = ya
    opti.subject_to(xa(1,j) >=ymin - delta_ipo(j));
    opti.subject_to(xa(1,j) <=ymax + delta_iper(j));
    %x2 e x3
    opti.subject_to(o3*(xa(2,j) + xa(3,j)) <= IOB_bound_eq(j) + delta_IOB_eq(j));
    %u
    opti.subject_to(ua(j) <=umax);
    opti.subject_to(ua(j) >=umin);
    %delta ipo e iper >=0
    opti.subject_to(delta_ipo(j) >=0);
    opti.subject_to(delta_iper(j) >=0);
    opti.subject_to(delta_IOB_eq(j) >=0);

    %opti.subject_to(xa(:,j) >=0);
end

opti.subject_to(xa(1,T+1) >=ymin - delta_ipo(T+1));
opti.subject_to(xa(1,T+1) <=ymax + delta_iper(T+1));
%x2 e x3
opti.subject_to(o3*(xa(2,T+1) + xa(3,T+1)) <= IOB_bound_eq(T+1) + delta_IOB_eq(T+1));
%delta ipo e iper >=0
opti.subject_to(delta_ipo(T+1) >=0);
opti.subject_to(delta_iper(T+1) >=0);
opti.subject_to(delta_IOB_eq(T+1) >=0);
%opti.subject_to(xa(:,T+1) >=0);


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
MPC=opti.to_function('MPC',{xk,rk,IOB_bound,Si_Tar_p,G_basal_p,IOB_basal_p,Si_Tar_eq,G_basal_eq,IOB_bound_eq},{uf,xf,V,xa,ua,delta_ipo,delta_iper,delta_ipo_x1,delta_iper_x1,delta_IOB,delta_IOB_eq,delta_XN,delta_simm},{'xk','rk','IOB_bound','Si_Tar_p','G_basal_p','IOB_basal_p','Si_Tar_eq','G_basal_eq','IOB_bound_eq'},{'uf_opt','xf_opt','V_opt','xa','ua','delta_ipo','delta_iper','delta_ipo_x1','delta_iper_x1','delta_IOB','delta_IOB_eq','delta_XN','delta_simm'});


end



