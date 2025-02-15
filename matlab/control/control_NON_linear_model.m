function control_NON_linear_model(data_circadian,patient,theta_ott_NL,Ts,deltaT)
%% IMPORT CASADI
import casadi.*

%% ESTRAZIONE DATI
[~,Gb,CR,CF,Ub,~,~,~,~,~] = data_extraction(data_circadian,patient);

%% DEFINIZIONE SISTEMA LINEARE
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

%% COSTANTI MPC
%orizzionte predittivo
N=72;
T=(60*24)/Ts;

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
ymin = 70;
ymax = 180;

%Pesi ipoglicemia e iperglicemia
p_ipo = 1*10^9;
p_iper = 1*10^7;
%Pesi variabili slack stati x1 e IOB(x2 e x3)
p_ipo_x1 =1*10^9;
p_iper_x1 =1*10^9;
p_IOB = 1*10^9;

p_XN = eye(4)*1*10^9;
p_SIM=1*10^9;

%% VARIABILI e PARAMETRI OTTIMIZZAZIONE
opti = casadi.Opti();

%variabili di ottimizzazione
xf = opti.variable(6,N+1);
uf = opti.variable(N);

%var. ausiliarie
xa = opti.variable(4,T+1);
ua = opti.variable(1,T);

%costo addizionale variabili slack
delta_iper_x1 = opti.variable(N);
delta_ipo_x1 = opti.variable(N);
delta_IOB = opti.variable(N);

%costo addizionale
delta_ipo = opti.variable(T);
delta_iper = opti.variable(T);
delta_IOB_eq = opti.variable(T);

delta_XN = opti.variable(4,1);
delta_simm = opti.variable();

%vincolo dinamico IOB(t)
IOB_bound = opti.parameter(N);
xk = opti.parameter(6,1); % x0 stimato da osservatore
rk = opti.parameter(N);

%parte NON LINEARE
%parametri
Si_Tar_p = opti.parameter(N);
G_basal_p = opti.parameter(N);
IOB_basal_p = opti.parameter();
%parametri per traiettoria equilibrio
Si_Tar_eq = opti.parameter(T);
G_basal_eq = opti.parameter(T);
IOB_bound_eq = opti.parameter(T);

%% FUNZIONE DI COSTO
Q=1;
R=100;
S=10^3;

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

V_XN = delta_XN(1:4)'*p_XN*delta_XN(1:4);

V_SIM = delta_simm^2*p_SIM;

V = V + V_XN + V_SIM;

opti.minimize(V);

%warm start ua
opti.set_initial(ua, Ub * ones(1, T));

%% VINCOLI

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
end

%intersezione traiettoria sistema e traiettoria equilibrio
% opti.subject_to( abs(xf(1:3,N+1) - xa(1:3,N+1)) <= delta_XN(1:3));
% opti.subject_to( abs(xf(6,N+1) - xa(4,N+1)) <= delta_XN(4));
opti.subject_to(xf(1:3,N+1) == xa(1:3,N+1));
opti.subject_to(xf(6,N+1) == xa(4,N+1));
%stato iniziale e terminale traiettoria eq devono coincidere
opti.subject_to(xa(:,1) == xa(:,T+1));

% opti.subject_to(xa(1:3,1) == xk(1:3));
% opti.subject_to(xa(4,1) == xk(6));
% opti.subject_to(abs(xa(:,1) - xa(:,T+1))<= delta_simm);

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
end

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
MPC=opti.to_function('MPC',{xk,rk,IOB_bound,Si_Tar_p,G_basal_p,IOB_basal_p,Si_Tar_eq,G_basal_eq,IOB_bound_eq},{uf,xf,V,xa,ua,delta_ipo,delta_iper,delta_ipo_x1,delta_iper_x1,delta_IOB,delta_IOB_eq,delta_XN},{'xk','rk','IOB_bound','Si_Tar_p','G_basal_p','IOB_basal_p','Si_Tar_eq','G_basal_eq','IOB_bound_eq'},{'uf_opt','xf_opt','V_opt','xa','ua','delta_ipo','delta_iper','delta_ipo_x1','delta_iper_x1','delta_IOB','delta_IOB_eq','delta_XN'});

%% OSSERVATORE
%KF settings
ki=1*10^-2;
kd=1;
%varianza
sigma_q = 1;
sigma_q_cgm = 6.51;
%Jacobiana sistema non lineare passo 0
A_j = [-o1 -CF 0 o2 0 -Ub;
    0 -1/o3 1/o3 0 0 0;
    0 0 -1/o3 0 0 0;
    0 0 0 -1/o4 1/o4 0;
    0 0 0 0 -1/o4 0;
    -o1/o6 -o3/o7 -o3/o7 0 0 -1/o5]; %CF=x0(6) Ub=x0(2)
%Discretizzazione jacobiana sistema non lineare
A_j_d = eye(size(A_j,1)) + Ts*A_j;

%Ingressi (matrice jacobiana rispetto all'ingresso u,r,Si)
B_u = [ 0 0 1/o3 0 0 0]';
B_r = [0 0 0 0 1/o4 0]';
%considero x6=Si un ingresso
B_si = [0 0 0 0 0 1/o5]';

E= [o0 0 0 0 0 0]';
C_d=[1 0 0 0 0 0];

%Discretizzazione matrici ingressi
B_u_d= Ts*B_u;
B_r_d = Ts*B_r;
B_si_d = Ts*B_si;
E_d = Ts*E;

Bd=[B_u_d B_r_d E_d B_si_d];

B_0 = [1 0 0 0 0 0]';

% Matrici sistema aumentato
Aaug = [A_j_d B_0 zeros(size(A_j_d,1),1);
    zeros(1,size(A_j_d,1)), 1, Ts;
    zeros(1,size(A_j_d,1)), 0, 1];

Baug = [Bd;zeros(2,size(Bd,2))];
G = [0 ki ki 0 0 0 0 kd]';
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
v_delta_iper_x1 = [];
v_delta_ipo_x1 = [];
v_delta_IOB =[];
v_delta_IOB_eq =[];
v_delta_XN =[];
v_xk_obs =[];
v_x_real =[];
v_time_sol =[];

%tempo di simulazione in minuti
Tmax = 24*60*1;

%ipotesi iniziali
x0 = [Gb Ub Ub 0 0 CF]';
u0 = 0;

%costruzione upper-lower bound IOB
IOB_vet = create_IOB_vector_NL(Tmax+(T*Ts),IOB_s,IOB_d);

%simulazione pasti (1 giorno + 6h)
% rk_in = zeros(1,Tmax+(N*Ts));
[rk_in,~] = create_RK_random(Tmax+(N*Ts));
% rk_in(481:510) = 60/30; %60g divisi in 30 minuti
% rk_in(111:140) = 60/30; %60g divisi in 30 minuti
% rk_in(121:150) = 60/30; %60g divisi in 30 minuti

%insulina da iniettare
u_in = 0;
%intervallo simulazione ode45
tspan = [0 1];

%Inizializzazioni parte NON lineare
Si_tar = CF*(1+0.01*o8*sin((2*pi*(1:(Tmax+T*Ts))*deltaT)/(60*24) + 2*pi*0.01*o9));
G_basal = (o0-Si_tar*Ub)/o1;
IOB_basal=2*Ub*o3;

%inizializzazioni osservatore
x_aug = [x0;0;0];
u_aug = [u0;0;1;Si_tar(1)];
y=Gb;
%inizializzazione primo stato sistema reale
x_real = x0;

for k=1:Tmax
    if mod(k-1,Ts) == 0
        %OSSERVATORE
        %stima x0 con osservatore
        [xk_obs,Pkk,Kk]=ODO(Aaug,Baug,Caug,y,x_aug,u_aug,Pkk,Q_kf,R_kf);
        %non considero i disturbi per MPC (prendo solo i primi 6 stati)
        xk_sim = xk_obs(1:6);

        %costruzione vettore rk_sim di lunghezza N (predizione 6h)
        rk_sim = create_RK_Sim(rk_in,k,Ts,N);

        %MPC
        %risolvo il problema di minimizzazione al passo k
        [uf_sol,xf_sol,V_sol,xa_sol,ua_sol,d_ipo,d_iper,d_ipo_x1,d_iper_x1,d_IOB,d_IOB_eq,d_XN]=MPC(xk_sim,rk_sim,IOB_vet(k:Ts:k+N*Ts-1),Si_tar(k:Ts:k+N*Ts-1),G_basal(k:Ts:k+N*Ts-1),IOB_basal,Si_tar(k:Ts:k+T*Ts-1),G_basal(k:Ts:k+T*Ts-1),IOB_vet(k:Ts:k+T*Ts-1));
        % [uf_sol,xf_sol,V_sol,d_ipo,d_iper,xa_sol]=MPC(xk_sim,rk_sim,IOB_vet(k:Ts:k+N*Ts-1),Si_tar(k:Ts:k+N*Ts-1),G_basal(k:Ts:k+N*Ts-1),IOB_basal,G_basal(k:Ts:k+T*Ts-1),Si_tar(k:Ts:k+T*Ts-1),IOB_vet(k:Ts:k+T*Ts-1));

        %converto risultato simbolico in risultato numerico
        xf_sim=full(xf_sol);
        uf_sim=full(uf_sol);
        Vn=full(V_sol);
        xa_sim = full(xa_sol);
        ua_sim = full(ua_sol);
        delta_iper_sim = full(d_iper);
        delta_ipo_sim = full(d_ipo);
        delta_iper_x1_sim = full(d_iper_x1);
        delta_ipo_x1_sim = full(d_ipo_x1);
        delta_IOB_sim = full(d_IOB);
        delta_IOB_eq_sim = full(d_IOB_eq);
        delta_XN_sim = full(d_XN);
        %verifico exit MPC
        exit=get_stats(MPC);
        if exit.success==1
            % prendo solo il primo elemento di uf ottima e lo applico al sistema
            u_in=uf_sim(1);
        else
            %stampo grafico soluzione precedente
            %TEST
            figure('Name', ['Soluzione precedente X']);
            plot(v_xa(end-3,:),'DisplayName', 'xa');
            hold on;
            legend('show');
            plot(v_xf(end-5,:),'DisplayName', 'xf');
            plot(G_basal(k:Ts:k+T*Ts-1),'DisplayName', 'Gb(t)');

            figure('Name', ['Soluzione precedente IOB']);
            plot(o3*(v_xa(end-2,:) + v_xa(end-1,:)) ,'DisplayName', 'xa');
            hold on;
            legend('show');
            plot(o3*(v_xf(end-4,:) + v_xf(end-3,:)),'DisplayName', 'xf');
            plot(IOB_vet(k:Ts:k+T*Ts-1),'DisplayName', 'IOB(t)');

            figure('Name', ['Soluzione precedente U']);
            stem(v_ua(end,:),'Marker','none','DisplayName', 'ua');
            hold on;
            legend('show');
            stem(v_uf(end,:),'Marker','none','DisplayName', 'uf');

            %stampo grafico soluzione ATTUALE NON FATTIBILE
            %TEST
            figure('Name', ['Soluzione attuale NON fattibile X']);
            plot(xa_sim(1,:),'DisplayName', 'xa');
            hold on;
            legend('show');
            plot(xf_sim(1,:),'DisplayName', 'xf');
            plot(G_basal(k:Ts:k+T*Ts-1),'DisplayName', 'Gb(t)');

            figure('Name', ['Soluzione attuale NON fattibile IOB']);
            plot(o3*(xa_sim(2,:) + xa_sim(3,:)) ,'DisplayName', 'xa');
            hold on;
            legend('show');
            plot(o3*(xf_sim(2,:) + xf_sim(3,:)),'DisplayName', 'xf');
            plot(IOB_vet(k:Ts:k+T*Ts-1),'DisplayName', 'IOB(t)');

            figure('Name', ['Soluzione attuale NON fattibile U']);
            stem(ua_sim(:),'Marker','none','DisplayName', 'ua');
            hold on;
            legend('show');
            stem(uf_sim(:),'Marker','none','DisplayName', 'uf');

            %se non trovo una soluzione ottima ingresso 0
            u_in=0;
        end

        %matrice stati x controllati con MPC su N=72 (6 ore)
        v_xf = [v_xf;xf_sim];
        %matrice stati x traiettoria equilibrio con MPC
        v_xa = [v_xa;xa_sim];
        %matrice controllo u MPC su N=72 (6 ore)
        v_uf = [v_uf;uf_sim'];
        %matrice controllo ua traiettoria equilibrio con MPC
        v_ua = [v_ua;ua_sim];
        %matrice valore del costo
        v_VN = [v_VN; Vn];
        %matrice soluzione MPC trovata
        v_exit=[v_exit,exit.success];
        %salvo matrici osservatore per vedere se convergono
        v_kk = [v_kk;Kk];
        v_pkk = [v_pkk;Pkk];
        %salvo valori disturbi per vedere se convergono
        v_d1 = [v_d1;xk_obs(7)];
        v_d2 = [v_d2;xk_obs(8)];
        %salvo valori delta calcolati in funzione MPC
        v_delta_iper = [v_delta_iper;delta_iper_sim];
        v_delta_ipo = [v_delta_ipo;delta_ipo_sim];
        v_delta_iper_x1 = [v_delta_iper_x1;delta_iper_x1_sim];
        v_delta_ipo_x1 = [v_delta_ipo_x1;delta_ipo_x1_sim];
        v_delta_IOB = [v_delta_IOB;delta_IOB_sim];
        v_delta_IOB_eq = [v_delta_IOB_eq;delta_IOB_eq_sim];
        v_delta_XN = [v_delta_XN;delta_XN_sim];
        %salvo valori xk osservati dall'osservatore
        v_xk_obs = [v_xk_obs xk_obs];
        %salvo tempi risoluzione MPC di ogni iterazione
        v_time_sol = [v_time_sol;exit.t_wall_total];

        %aggiorno valori aug per prossima iterazione MPC
        x_aug = xk_obs;
        u_aug = [u_in;rk_in(k);1;Si_tar(k)];
        %aggiorno valore matrice Aaug per iterazione sucessiva
        Aaug = updateAaug(xk_obs,A_j,B_0,Ts);
    end

    %avanzamento sistema reale
    [tt,x_sim] = ode45(@(t, x) patient_fun_NL(t, x, theta_ott_NL, G_basal(k),IOB_basal,Si_tar(k), u_in, rk_in(k)), tspan, x_real);
    % Prendo l'ultimo stato
    x_real = x_sim(end, :)';
    %unica misura che sono in grado di effettuare è x1
    y=x_real(1);

    %salvo i valori insulina inserita
    v_u= [v_u,u_in];

    %azzero ingresso insulina per ciclo sucessivo
    %u_in = 0;
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
save('dati_intermedi_NL.mat', 'v_xf', 'v_xa','v_ua','v_uf','v_u', 'v_y', 'v_exit', 'v_VN', ...
    'v_x_sim', 'v_tt', 'v_pkk', 'v_kk', 'v_d1', 'v_d2', ...
    'v_delta_ipo', 'v_delta_iper', 'v_xk_obs','v_x_real','v_time_sol','IOB_vet','rk_in');

%% GRAFICI
figure('Name', ['Controllo MPC - Paziente ' num2str(patient)]);

%---------GLICEMIA---------
subplot(5, 1, 1);
plot(v_y(1:Ts:end), 'r-', 'LineWidth', 1, 'DisplayName', 'Glicemia reale');
hold on;
plot(v_xk_obs(1,:), 'b-', 'LineWidth', 1, 'DisplayName', 'Glicemia osservatore (ODO)');
plot(G_basal(1:Ts:end), 'm--', 'LineWidth', 1, 'DisplayName', 'Gbasal(t)');
hold off;
grid on;
xlim([1, length(v_xk_obs)]);
ylabel('Glicemia [mg/dL]');
legend('show');
set(gca, 'FontSize', 12);

%------------IOB------------
IOB_obs = o3*(v_xk_obs(2,:) + v_xk_obs(3,:));
IOB_real = o3*(v_x_real(2,:) + v_x_real(3,:));
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
Ra_real = o2*v_x_real(4,:);
Ra_obs = o2*v_xk_obs(4,:);
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

clc;
disp(['tempo totale: ' num2str(sum(v_time_sol)/60)]);