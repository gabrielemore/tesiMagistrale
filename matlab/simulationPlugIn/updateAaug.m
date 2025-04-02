function Aaug = updateAaug(xk_obs,A_j,B_0,Ts)
%update x6 e x2 (il resto sono parametri costanti)
A_j(1,2) = -xk_obs(6);
A_j(1,6) = -xk_obs(2);

%discretizzo
A_j_d = eye(size(A_j,1)) + Ts*A_j;

Aaug = [A_j_d B_0 zeros(size(A_j_d,1),1);
    zeros(1,size(A_j_d,1)), 1, Ts;
    zeros(1,size(A_j_d,1)), 0, 1];
end
