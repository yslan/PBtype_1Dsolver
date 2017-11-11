function res = TotalRes_Poisson(Phi,C1_ini,C2_ini)
% Given Phi, compute total res of PDE in form of Phi

global NN TotNumDM DegDM
global Xprolong xmin xmax
global MassVec JacVec Diff
global k_B T e zval
global Psi_1 Psi_2


[C1,C2] = solve_c_by_phi(Phi,C1_ini,C2_ini,'fsolve');

% Boundary condition-> force to be Dirichlet 
pen_BC.alpha_m = 1; pen_BC.beta_m = 0;  % specifies here
pen_BC.alpha_p = 1; pen_BC.beta_p = 0;  % specifies here


Bdry_D = [-1;0];
Bdry_N = [0;0];


vec_one = ones(size(Xprolong));

% construct ML's
    [~,DegDM,BlockL,SupMat,SubMat,Fprolong] =...
                get_L_Fvec(xmin,xmax,NN,TotNumDM,pen_BC,Bdry_D,Bdry_N,...
                           vec_one,0*vec_one);


wk_rhs = MassVec .* (zval(1)*C1 + zval(2)*C2);


ML_Phi = MLx(Phi,DegDM,BlockL,SupMat,SubMat) - Fprolong;

res = ML_Phi./JacVec - wk_rhs;

end
