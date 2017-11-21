function [f1,f2] = userf(C1,C2)
global param

global TotNumDM NN xmin xmax
global k_B T e zval epsilon
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff

% f1 = zeros(size(C1));
% f2 = zeros(size(C2));

rad_1 = 1; % unit: 1E-10 m 
rad_2 = 1;

rad3_1= rad_1 ^3;
rad3_2= rad_2 ^3;

f1 = k_B*T*log(rad3_1) - k_B*T*log(1-rad3_1*C1-rad3_2*C2);
f2 = k_B*T*log(rad3_2) - k_B*T*log(1-rad3_1*C1-rad3_2*C2);


end
