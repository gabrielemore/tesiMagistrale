function IOB_vet = create_IOB_vector(minuti,IOB_s,IOB_d)
%orizzonte temporale totale +6h predizione
m_tot = minuti+360;
%inizializzo IOB_vet a 0
IOB_vet = zeros(1,m_tot);

for i=1:m_tot
    %se orario giorno (6h-22h] -> IOB_d
    if mod(i,1440) > 360 && mod(i,1440)<=1320
        IOB_vet(i)=IOB_d;
        %altrimenti ->IOB_s
    else
        IOB_vet(i)=IOB_s;
    end
end


% IOB_vet = zeros(1,1800);
% IOB_vet(1:360) = IOB_s;
% IOB_vet(361:1320) = IOB_d;
% IOB_vet(1321:1440) = IOB_s;
% % + 6h
% IOB_vet(1441:1800) = IOB_s;





