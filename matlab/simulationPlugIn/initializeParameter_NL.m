function [o0,o1,o2,o3,o4,o5,o6,o7,o8,o9] = initializeParameter_NL(patient_name)

patient_number = sscanf(patient_name,'adult#%d');
file_id = sprintf('p%03d', patient_number);
file_name = sprintf('PARAMETRI_STIMATI/NL_CIRC_ON/theta_ott_non_lineare_%s.mat', file_id);

load(file_name);

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
end

