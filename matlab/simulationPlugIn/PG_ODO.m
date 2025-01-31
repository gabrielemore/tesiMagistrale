function [x_aug,Pkk,Kk] = PG_ODO(Aaug,Baug,Caug,y,x_aug,u_aug,Pkk,Q_kf,R_kf)

%P(k|k-1)
Pk = Aaug*Pkk*Aaug' + Q_kf;
%K(k)
Kk = Pk*Caug'*inv(Caug*Pk*Caug'+ R_kf);
%x_cap(k|k)
x_aug = Aaug*x_aug + Baug*u_aug + Kk *(y - Caug*(Aaug*x_aug + Baug*u_aug));
%P(k|k) e al ciclo sucessivo sarà P(k-1|k-1)
Pkk = (eye(size(Kk,1)) - Kk*Caug)*Pk;

end

