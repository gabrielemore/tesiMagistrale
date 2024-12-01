function [J] = cost_function_non_linear_system(theta,x0,y,u,r,Ts,T,theta0)

o0=theta(1);
o1=theta(2);
o2=theta(3);
o3=theta(4);
o4=theta(5);
o5=theta(6);
o6=theta(7);
o7=theta(8);
o8=theta(9);
o9=theta(10);

%IOB basal
IOB_basal=2*x0(2)*o3; %x0(2)=Ub

%theta^0 per funzione di minimizzazione
o_0 = [o2-theta0(3), o3-theta0(4), o4-theta0(5)];
%alfa per funzione di minimizzazione
alfa = diag([1,1,1/theta0(5)^2]);
%regularizaiton term
reg= o_0*alfa*o_0';
   
x(:,1) = x0;
y_cap = zeros(T,1);
deltaG = zeros(T,1);
deltaIOB = zeros(T,1);
IOB = zeros(T,1);

i = 1:T;
%baseline insulin sensitivity
Si_tar = x0(6)*(1+0.01*o8*sin((2*pi*i*Ts)/(60*24) + 2*pi*0.01*o9)); % x0(6)=CF
%non abbiamo G_b fisso ma variabile 
G_basal = (o0-Si_tar*x0(2))/o1; %x0(2)=Ub;

for k=1:T
    %deviation of gluscose from its basal
    deltaG(k) = x(1,k) - G_basal(k);

    %deviation of IOB from its basal
    IOB(k) = o3*(x(2,k) + x(3,k));
    deltaIOB(k) = IOB(k) - IOB_basal;

    %sistema non lineare
    x_dot(1)= o0 - (o1*x(1,k)) - (x(6,k)*x(2,k)) + (o2*x(4,k)); %G
    x_dot(2)= -(1/o3 * x(2,k)) + (1/o3 * x(3,k)); %Qi
    x_dot(3)= -(1/o3 * x(3,k)) + (1/o3 * u(k)); %Qisub
    x_dot(4)= -(1/o4 * x(4,k)) + (1/o4 * x(5,k)); %Qg
    x_dot(5)= -(1/o4 * x(5,k)) + (1/o4 * r(k)); %Qsto
    x_dot(6)= -(1/o5 * x(6,k)) -(o1/o6 * deltaG(k)) - (1/o7 * deltaIOB(k)) + (1/o5 * Si_tar(k)); %Si

    %discretizzazione
    x(:,k+1) = x(:,k) + Ts*x_dot(:);
    
    y_cap(k) = x(1,k);
   
end

J = (norm(y(1:T,1)-y_cap))^2/(norm(y(1:T,1)-mean(y(1:T))))^2 + reg;
