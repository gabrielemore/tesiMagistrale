function [N,umin,umax,x1min,x1max,CHO_UB,tau,IOB_s,IOB_d,ymin,ymax,p_ipo,p_iper,p_ipo_x1,p_iper_x1,p_IOB]=initializeConstantMPC(o4,Ub,CR)

%orizzionte predittivo
N=72;
%vincoli controllo
umin = 0;
umax = 15;

%vincoli (nb:x4 e x5 senza vincoli)
x1min = 54;
x1max = 300;
% IOB
CHO_UB = 90;
tau = 120;
IOB_s = o4*2*Ub; %[22 - 6)
IOB_d = IOB_s + (CHO_UB/CR + tau*Ub); %[6 - 22)

%uscita
ymin = 70;
ymax = 140;

%Pesi ipoglicemia e iperglicemia
p_ipo = 1*10^9;
p_iper = 1*10^7;
%Pesi variabili slack stati x1 e IOB(x2 e x3)
p_ipo_x1 =1*10^9;
p_iper_x1 =1*10^9;
p_IOB = 1*10^9;

end

