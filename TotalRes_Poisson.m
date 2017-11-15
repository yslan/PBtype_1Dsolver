function res = TotalRes_Poisson(Phi,C1_ini,C2_ini)
% Given Phi, compute total res of PDE in form of Phi

global param

global TotNumDM NN xmin xmax
global k_B T e zval epsilon
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff

[C1,C2] = solve_c_by_phi(Phi,C1_ini,C2_ini,'fsolve');

% Boundary condition-> force to be Dirichlet 
pen_BC.alpha_m = 1; pen_BC.beta_m = 0;  % specifies here
pen_BC.alpha_p = 1; pen_BC.beta_p = 0;  % specifies here


Bdry_D = [phi_L;0];
Bdry_N = [0;0];


vec_one = ones(size(Xprolong));

% construct ML's
    [~,DegDM,BlockL,SupMat,SubMat,Fprolong] =...
                get_L_Fvec(xmin,xmax,NN,TotNumDM,pen_BC,Bdry_D,Bdry_N,...
                           epsilon*vec_one,0*vec_one);


wk_rhs = MassVec .* (zval(1)*C1 + zval(2)*C2);


ML_Phi = MLx(Phi,DegDM,BlockL,SupMat,SubMat) - Fprolong;

res = ML_Phi./JacVec - wk_rhs;

end
