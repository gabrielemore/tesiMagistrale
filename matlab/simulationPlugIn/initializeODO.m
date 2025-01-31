function [Aaug,Baug,Caug,Q_kf,R_kf,Pkk] = initializeODO(A_d,B_u_d,B_r_d,C_d,E_d,Ts)
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

[~,Pkk,~,~] = dlqe(Aaug,G,Caug,sigma_q,R_kf);
end

