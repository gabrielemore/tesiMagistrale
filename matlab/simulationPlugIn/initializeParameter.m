function [o0,o1,o2,o3,o4,o5] = initializeParameter()
load('theta_ott_lineare.mat');

o0=theta_ott_ML(1);
o1=theta_ott_ML(2);
o2=theta_ott_ML(3);
o3=theta_ott_ML(4);
o4=theta_ott_ML(5);
o5=theta_ott_ML(6);
end

