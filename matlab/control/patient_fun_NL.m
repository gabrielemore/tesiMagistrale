function x_dot = patient_fun_NL(t,x,theta_ott_NL,G_basal,IOB_basal,Si_tar,u,r)

o0=theta_ott_NL(1);
o1=theta_ott_NL(2);
o2=theta_ott_NL(3);
o3=theta_ott_NL(4);
o4=theta_ott_NL(5);
o5=theta_ott_NL(6);
o6=theta_ott_NL(7);
o7=theta_ott_NL(8);
o8=theta_ott_NL(9);
o9=theta_ott_NL(10);

x_dot = zeros(6,1);

%sistema non lineare
x_dot(1)= o0 - (o1*x(1)) - (x(6)*x(2)) + (o2*x(4)); %G
x_dot(2)= -(1/o3 * x(2)) + (1/o3 * x(3)); %Qi
x_dot(3)= -(1/o3 * x(3)) + (1/o3 * u); %Qisub
x_dot(4)= -(1/o4 * x(4)) + (1/o4 * x(5)); %Qg
x_dot(5)= -(1/o4 * x(5)) + (1/o4 * r); %Qsto
x_dot(6)= -(1/o5 * x(6)) -(o1/o6 * (x(1) - G_basal)) - (1/o7 * (o3*(x(2) + x(3)) - IOB_basal)) + (1/o5 * Si_tar); %Si
end

