function IOB_vet = PGcreate_IOB_vector(minuti,IOB_s,IOB_d)
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






