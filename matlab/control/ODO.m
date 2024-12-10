function [x_aug,Pkk] = ODO(Aaug,Baug,Caug,y,x_aug,u_aug,Pkk,Q_kf,R_kf)

%P(k|k-1)
Pk = Aaug*Pkk*Aaug' + Q_kf;
%K(k)
Kk = Pk*Caug'*(Caug*Pk*Caug'+ R_kf)^-1;
%x_cap(k|k)
x_aug = Aaug*x_aug + Baug*u_aug + Kk *(y - Caug*x_aug);
%P(k|k) e al ciclo sucessivo sarebbe P(k-1|k-1)
Pkk = (eye(7) - Kk*Caug)*Pk;

end

