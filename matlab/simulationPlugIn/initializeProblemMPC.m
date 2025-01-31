function [opti,MPC] = initializeProblemMPC(N,A_d,B_u_d,B_r_d,C_d,E_d,x1min,x1max,p_IOB,p_iper_x1,p_ipo_x1,p_iper,p_ipo,umax,umin,ymin,ymax,o4)

opti = casadi.Opti();

%variabili di ottimizzazione
xf = opti.variable(5,N+1);
uf = opti.variable(N);

%var. ausiliarie
xa = opti.variable(3,1);
ua = opti.variable();

%costo addizionale variabili slack
delta_iper_x1 = opti.variable(N);
delta_ipo_x1 = opti.variable(N);
delta_IOB = opti.variable(N);

%costo addizionale
delta_ipo = opti.variable();
delta_iper = opti.variable();

%vincolo dinamico IOB(t)
IOB_bound = opti.parameter(N+1);
xk = opti.parameter(5,1); % x0 stimato da osservatore
rk = opti.parameter(N);

Q=1;
R=100;

V=0;
for j=1:N
    %funzione costo
    V_dyn = (C_d*xf(:,j) - xa(1))'*Q*(C_d*xf(:,j) - xa(1)) + (uf(j) - ua)'*R*(uf(j) - ua);
    %costo slack IOB
    V_IOB = p_IOB*delta_IOB(j)^2;
    %costo slack x1
    V_x1 = p_iper_x1*delta_iper_x1(j)^2 + p_ipo_x1*delta_ipo_x1(j)^2;
    %costo tot
    V = V + V_dyn + V_IOB + V_x1;
end
%costo addizionale
V_s = p_iper*delta_iper^2 + p_ipo*delta_ipo^2;
%costo totale
V = V + V_s;

opti.minimize(V);

%condizione iniziale x0
opti.subject_to(xf(:,1)==xk);

for j=1:N
    %dinamica modello
    opti.subject_to(xf(:,j+1)==A_d*xf(:,j) + B_u_d*uf(j) + B_r_d*rk(j) + E_d);
    %x1 +- slack
    opti.subject_to(xf(1,j) >=x1min - delta_ipo_x1(j));
    opti.subject_to(xf(1,j) <=x1max + delta_iper_x1(j));
    %x2 e x3 +- slack
    %opti.subject_to(o4*(xf(2,j)+xf(3,j)) >=0);
    opti.subject_to(o4*(xf(2,j)+xf(3,j)) <=IOB_bound(j) + delta_IOB(j));
    %u
    opti.subject_to(uf(j) <=umax);
    opti.subject_to(uf(j) >=umin);
    %slack x1 e IOB >=0
    opti.subject_to(delta_ipo_x1(j) >=0);
    opti.subject_to(delta_iper_x1(j) >=0);
    opti.subject_to(delta_IOB(j) >=0);
end

%eq
%x1a = ya
opti.subject_to(xa(1) >=ymin - delta_ipo);
opti.subject_to(xa(1) <=ymax + delta_iper);
%delta ipo e iper >=0
opti.subject_to(delta_ipo >=0);
opti.subject_to(delta_iper >=0);
%x2 e x3
%opti.subject_to(o4*(xa(2)+xa(3)) >=0);
opti.subject_to(o4*(xa(2) + xa(3)) <= IOB_bound(N+1));
%ua
opti.subject_to(ua <=umax);
opti.subject_to(ua >=umin);
%C_n*x(N)=xa
opti.subject_to(xf(1:3,N+1)==xa);
opti.subject_to(xa==A_d(1:3,1:3)*xa + B_u_d(1:3)*ua + E_d(1:3));

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
MPC=opti.to_function('MPC',{xk,rk,IOB_bound},{uf,xf,V,delta_ipo,delta_iper},{'xk','rk','IOB_bound'},{'uf_opt','xf_opt','V_opt','delta_ipo','delta_iper'});

end

