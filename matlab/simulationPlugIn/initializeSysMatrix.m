function [A_d,B_u_d,B_r_d,C_d,E_d] = initializeSysMatrix(o0,o1,o2,o3,o4,o5,Ts,deltaT)

%Forma matriciale
A = [-o1 -o2 0 o3 0;0 -1/o4 1/o4 0 0;0 0 -1/o4 0 0; 0 0 0 -1/o5 1/o5; 0 0 0 0 -1/o5];
B_u= [0 0 1/o4 0 0]';
B_r= [0 0 0 0 1/o5]';
E=[o0 0 0 0 0]';
C=[1 0 0 0 0];
%discretizzazione
[A_d,B_u_d,B_r_d,C_d,E_d]=PGsystem_discretization(A,B_u,B_r,E,C,Ts,deltaT);
end

