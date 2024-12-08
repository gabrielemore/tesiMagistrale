clear all
close all
clc

warning off all

%% importo le librerie di Casadi
addpath ('/Users/af/Documents/MATLAB/casadi-osx-v3.6.3')
%addpath ('/Users/af/Documents/MATLAB/casadi-M1-v3.6.3')
import casadi.*

%% Model

A=[1 1; 0 1];
B=[0 0.5; 1 0.5];
C=[1 0];
D=[0 0];

[nx,nu]=size(B);
ny=length(C(:,1));

xmax=[5;5]; % vincoli di massimo e minimo sullo stato
xmin=-xmax;

umax=[.5; .5]; %vincoli di massimo e minimo sull'ingresso
umin=-umax;

%% dichiero stati e ingressi come variabili dell'intorno casadi.

x1=MX.sym('x1'); 
x2=MX.sym('x2');

x=[x1;x2];

u1=MX.sym('u1');
u2=MX.sym('u2');

u=[u1;u2];

%% Scrivo il modello in spazio di stati come funzione in Casadi
x_next=A*x+B*u;
F=Function('F',{x,u},{x_next},{'x','u'},{'x_next'}); %modelo discreto 

% % se avessi un modello tempo continuo scriverei questo
% ode=A*+B*u;
% 
% % continuous model
% f=Function('f',{x,u},{ode},{'x','p'},{'ode'});
% 
% % discretize model with rk
% Ts=1; % tempo di campionamento, da prendere dal paper
% intg_options=struct;
% intg_options.tf=Ts;
% intg_options.simplify=true;
% intg_options.number_of_finite_elements=5;
% 
% %DAE problem model
% dae=struct;
% dae.x=x;
% dae.p=u;
% dae.ode=f(x,u);
% intg=integrator('intg','rk',dae,intg_options);
% 
% res=intg('x0',x,'p',u); %simbolico
% x_next=res.xf;
% F=Function('F',{x,u},{x_next},{'x','p'},{'x_next'}); %modello discreto (un passo)

%% control problem
%opti = casadi.Opti('conic'); %necessaria per gurobi
opti = casadi.Opti();

N=15;

% optimzation variables
xf = opti.variable(nx,N+1);
uf = opti.variable(nu,N);

% xa e ua sono varibili di ottimizzazione ausiliari. Mi servono per ragioni
% di stabilità (come nel paper).

xa = opti.variable(nx,1); 
ua = opti.variable(nu,1);

ys = opti.parameter(ny,1);
xk = opti.parameter(nx,1);

%% Cost function

Q=eye(nx);
R=eye(nu);
T=100*eye(ny);

V=0;

for j=1:N
    
    L=(xf(:,j)-xa)'*Q*(xf(:,j)-xa)+(uf(:,j)-ua)'*R*(uf(:,j)-ua);
    V=V+L;
end

Vo=(C*xa-ys)'*T*(C*xa-ys); % costo di offset

V=V+Vo;
opti.minimize(V);

%% Constraints
opti.subject_to(xf(:,1)==xk);  %initial condition

for j=1:N
    opti.subject_to(xf(:,j+1)==F(xf(:,j),uf(:,j))); %predizione dinamica
    opti.subject_to(xf(:,j) >=xmin); 
    opti.subject_to(xf(:,j) <=xmax);
    opti.subject_to(uf(:,j) <=umax);
    opti.subject_to(uf(:,j) >=umin);
end

opti.subject_to(xa >=xmin);
opti.subject_to(xa <=xmax);
opti.subject_to(ua <=umax);
opti.subject_to(ua >=umin);
opti.subject_to(xf(:,N+1)==xa); % terminal constraint
opti.subject_to(xa==F(xa,ua));


%% solver

opti.solver('ipopt');
%opti.solver('gurobi');

% costruisco la funzione MPC per risolver il problema di ottimizzazione
MPC=opti.to_function('MPC',{xk,ys},{uf,xf,V,xa,ua},{'xk','ys'},{'uf_opt','xf_opt','V_opt','xa_opt','ua_opt'});

%% PROVA

vx=[];
vu=[];
vys=[];
vexit=[];
vVN=[];
Vn=0;


Tmax = 20; % tempo massimo di simulazione

% target
ys=0;

x0=[-2;0.3];
xk=x0;
uk=[0;0];

vys=[vys,ys];
vx=[vx,xk];

for k=1:Tmax
    uk_past=uk;

    [uf_sol,xf_sol,V_sol,xa_sol,ua_sol]=MPC(xk,ys);
   
    exit=get_stats(MPC);
    vexit=[vexit,exit.success];
  
%     if exit.success==1
%    
        uf=full(uf_sol); xf=full(xf_sol);
        Vn=full(V_sol);
%     else
%         uk=uk_past;
%     end
   
   % prendo solo il primo elemento di uf ottima e lo applico al sistema
   uk=uf(:,1);
   xk=full(F(xk,uk)); % faccio avanzare il sistema
   
   vys=[vys,ys];
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
