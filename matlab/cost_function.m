function J = cost_function(theta,x0,y,u,r,Ts,T,theta0)

o0=theta(1);
o1=theta(2);
o2=theta(3);
o3=theta(4);
o4=theta(5);
o5=theta(6);

CR=(theta0(4)^-1)*theta0(3);
%theta^0 per funzione di minimizzazione  !!!problema SEMPRE A 0
o_0 = [o2 - CR*o3;o4-theta0(5); o5-theta0(6)]';
%alfa per funzione di minimizzazione
alfa = diag([1,1,1/theta0(6)]);
%y average
y_m = mean(y);

%Definizione sistema lineare in forma matriciale
A = [-o1 -o2 0 o3 0;0 -1/o4 1/o4 0 0;0 0 -1/o4 0 0; 0 0 0 -1/o5 1/o5; 0 0 0 0 -1/o5];
B_u= [0 0 1/o4 0 0]';
B_r= [0 0 0 0 1/o5]';
E=[o0 0 0 0 0]';
C=[1 0 0 0 0];

J=0;
x(:,1)=x0;

for k=1:T
    x(:,k+1)=x(:,k) + Ts*[A*x(:,k) + B_u*u(k) + B_r*r(k) + E];
    y_cap(k)=C*x(:,k);
    J=J+ (norm(y(k)-y_cap(k))^2/norm(y(k)-y_m)^2 + o_0*alfa*o_0');
end

