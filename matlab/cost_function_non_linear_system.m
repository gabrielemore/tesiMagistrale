function J = cost_function_non_linear_system(theta,x0,y,u,r,Ts,T,theta0,CF)

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

%ottengo CR per trovare o_0 
CR_0=theta0(2)/theta0(3);
   
x_dot(:,1) = x0;
y_cap = zeros(T,1);

for k=1:T

    %sistema non lineare
    x_dot(1, k+1)= x_dot(1,k) + Ts*(o0 - o1*x_dot(1,k) - x_dot(6,k)*x_dot(2,k) + o2*x_dot(4,k)); %G
    x_dot(2, k+1)= x_dot(2,k) + Ts*(-1/o3 * x_dot(2,k+1) + 1/o3 * x_dot(3,k)); %Qi
    x_dot(3, k+1)= x_dot(3,k) + Ts*(-1/o3 * x_dot(3,k) + 1/o3 * u(k)); %Qisub
    x_dot(4, k+1)= x_dot(4,k) + Ts*(-1/o4 * x_dot(4,k) + 1/o4 * x_dot(5,k)); % Qg
    x_dot(5, k+1)= x_dot(5,k) + Ts*(-1/o4 * x_dot(5,k) + 1/o4 * r(k)); %Qsto
    
    %baseline insulin sensitivity
    Si_tar(k) = CF*(1+o8* sin((2*pi*k)/(60*24)) + 2*pi*o9);

    %deviation of gluscose from its basal
    %non abbiamo G_b fisso ma varia 
    G_basal(k) = (o0-Si_tar(k)*x0(2))/o1; %x0(2)=Ub;
    deltaG(k) = x_dot(1,k) - G_basal(k);

    %deviation of IOB from its basal
    IOB_basal=2*x0(2)*o3; %x0(2)=Ub
    IOB(k) = o3*(x_dot(2,k) + x_dot(3,k));
    deltaIOB(k) = IOB(k) - IOB_basal;
    
    x_dot(6, k+1)= x_dot(6,k) + Ts*(-1/o5 * x_dot(6,k) -1/o6 * deltaG(k) - 1/o7 * deltaIOB(k) + 1/o5 * Si_tar(k)); %Si
    
    y_cap(k) = x_dot(1,k);

    % costruzione vettore media al tempo k 
    y_m(k)=mean(y(1:k));

    %theta^0 per funzione di minimizzazione aggiornato ad ogni iterazione
    % con Si(k) al posto di o2 e gli altri elementi shiftati di 1
    o_0 = [x_dot(6,k)-CR_0*o2, o3-theta0(4), o4-theta0(5)];
    %alfa per funzione di minimizzazione
    alfa = diag([1,1,1/theta0(5)^2]);
    %regularizaiton term
    reg= o_0*alfa*o_0';

    if k==1
        J=0;
    else
        J = J + (norm(y(1:k)-y_cap(1:k)))^2/(norm(y(1:k)-mean(y(1:k))))^2 + reg;
    end
   
end

