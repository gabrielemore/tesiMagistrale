function rk_sim = create_RK_Sim(rk_in,k,Ts,N)
%Inizializzo il vettore simulazione a 0
rk_sim = zeros(1, N);
for i = 1:N
    %se ho un pasto al passo k
    if rk_in(k) ~= 0
        %Inserisco il pasto nel vettore della simulazione MPC fino a che
        %questo pasto non termina (incontro uno zero)
        rk_sim(i) = rk_in(k);
        %Incremento k con il passo temporale Ts
        k = k + Ts;
    else
        %Se ho uno zero (pasto terminato o mai iniziato) allora esco dal ciclo
        break;
    end
end
end