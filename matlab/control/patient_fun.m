function x_dot = patient_fun(t,x,A,B_u,B_r,E,u,r)
x_dot = zeros(5,1);
x_dot = A*x + B_u*u + B_r*r + E;
end

