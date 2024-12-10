function [A_d,B_u_d,B_r_d,C_d,E_d] = system_discretization(A,B_u,B_r,E,C,T,deltaT)

A_d = expm(A*T);
B_u_d = expm(A*(T-deltaT))*inv(A)*(expm(A*deltaT) - eye(5))*B_u;
B_r_d = inv(A)*(expm(A*T) - eye(5))*B_r;
C_d = C;
E_d = inv(A)*(expm(A*T) - eye(5))*E;

end

