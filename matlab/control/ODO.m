function [x_aug,Pkk] = ODO(Aaug,Baug,C_d,y,x_aug,u_aug,Pkk)

%P(k|k-1)
Pk = Aaug*Pkk*Aaug' + Q_kf;
%K(k)
Kk = Pk*C_d'*(C_d*Pk*C_d'+ R_kf)^-1;
%x_cap(k|k)
x_aug = Aaug*x_aug + Baug*u_aug + Kk *(y - C_d*x_aug);
%P(k|k) e al ciclo sucessivo sarebbe P(k-1|k-1)
Pkk = (eye(7) - Kk*C_d)*Pk;

end

