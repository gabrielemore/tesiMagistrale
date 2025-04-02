function [o0,o1,o2,o3,o4,o5] = initializeParameter(patient_name)

patient_number = sscanf(patient_name,'adult#%d');
file_id = sprintf('p%03d', patient_number);
%file_name = sprintf('PARAMETRI_STIMATI/LIN_CIRC_OFF/theta_ott_lineare_%s.mat', file_id);
file_name = sprintf('PARAMETRI_STIMATI/LIN_CIRC_ON/theta_ott_lineare_circ_on_%s.mat', file_id);

load(file_name);

% load('PARAMETRI_STIMATI/theta_ott_lineare_p001.mat');

% o0=theta_ott_ML(1);
% o1=theta_ott_ML(2);
% o2=theta_ott_ML(3);
% o3=theta_ott_ML(4);
% o4=theta_ott_ML(5);
% o5=theta_ott_ML(6);
o0=theta_ott_ML_circ_ON(1);
o1=theta_ott_ML_circ_ON(2);
o2=theta_ott_ML_circ_ON(3);
o3=theta_ott_ML_circ_ON(4);
o4=theta_ott_ML_circ_ON(5);
o5=theta_ott_ML_circ_ON(6);
end

