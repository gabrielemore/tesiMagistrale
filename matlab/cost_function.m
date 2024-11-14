function J = cost_function(theta,x0,y,u,r,Ts,T,theta0)

o0=theta(1);
o1=theta(2);
o2=theta(3);
o3=theta(4);
o4=theta(5);
o5=theta(6);

%ottengo CR per trovare o_0 
CR_0=theta0(3)/theta0(4);
%theta^0 per funzione di minimizzazione
o_0 = [o2-CR_0*o3, o4-theta0(5), o5-theta0(6)];
%alfa per funzione di minimizzazione
alfa = diag([1,1,1/theta0(6)^2]);
%regularizaiton term
reg= o_0*alfa*o_0';


%Definizione sistema lineare in forma matriciale
A = [-o1 -o2 0 o3 0;0 -1/o4 1/o4 0 0;0 0 -1/o4 0 0; 0 0 0 -1/o5 1/o5; 0 0 0 0 -1/o5];
B_u= [0 0 1/o4 0 0]';
B_r= [0 0 0 0 1/o5]';
E=[o0 0 0 0 0]';
C=[1 0 0 0 0];

%versione con theta2 e theta3 dimensionalmente sbagliate
%J=0;
%y_m = mean(y(1:T));

x(:,1)=x0;
y_cap = zeros(T,1);

for k=1:T
    %sistema lineare
    x(:,k+1) = x(:,k) + Ts*(A*x(:,k) + B_u*u(k) + B_r*r(k) + E);
    y_cap(k)=C*x(:,k);
    
    %versione con theta2 e theta3 dimensionalmente sbagliate
    % %y average al tempo k
    % y_m_k = mean(y(1:k));
    % J = J + (norm(y(k)-y_cap(k)))^2/(norm(y(k)-y_m))^2 + reg ;
    
    % costruzione vettore media al tempo k versione J fuori dal ciclo
    y_m(k)=mean(y(1:k));
    
    % %versione lenta, risultati simili a J fuori dal ciclo
    % if k==1
    %     J=0;
    % else
    %     J = J + (norm(y(1:k)-y_cap(1:k)))^2/(norm(y(1:k)-mean(y(1:k))))^2 + reg;
    % end
end

J = (norm(y(1:T,1)-y_cap))^2/(norm(y(1:T,1)-y_m'))^2 + reg;