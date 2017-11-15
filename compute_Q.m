function Q = compute_Q(C1,C2,mode)
% total charge
% mode: 0 = compute, 1 = exact
global param

global TotNumDM NN xmin xmax
global k_B T e zval epsilon
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff

z=  abs(zval(1)); % abs(z_1) = abs(z_2)

switch mode
    case 0
        Q = e * (MassVec.*JacVec)' * (zval(1)*C1 + zval(2)*C2);
    case 1
        Q = -2 * sqrt(2*k_B*T*epsilon*c0) * sinh(e*z*phi_L/(2*k_B*T));

end




end
