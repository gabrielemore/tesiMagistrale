function [Aaug,Baug,Caug,Q_kf,R_kf,Pkk,A_j,B_0] = initializeODO_NL(o0,o1,o2,o3,o4,o5,o6,o7,CF,Ub,Ts)
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
end

