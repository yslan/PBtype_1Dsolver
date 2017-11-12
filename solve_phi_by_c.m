function Phi = solve_phi_by_c(C1,C2,Phi_ini)
% Given C1 C2, solve Phi

global NN TotNumDM DegDM
global Xprolong xmin xmax
global MassVec JacVec Diff
global k_B T e zval
global c0 phi_L Psi_1 Psi_2

% Boundary condition-> force to be Dirichlet 
pen_BC.alpha_m = 1; pen_BC.beta_m = 0;  % specifies here
pen_BC.alpha_p = 1; pen_BC.beta_p = 0;  % specifies here


Bdry_D = [phi_L,0];
Bdry_N = [0,0];


vec_one = ones(size(Xprolong));
rhs = e*(zval(1)*C1 + zval(2)*C2);

% construct ML's
    [~,DegDM,BlockL,SupMat,SubMat,Fprolong] =...
                get_L_Fvec(xmin,xmax,NN,TotNumDM,pen_BC,Bdry_D,Bdry_N,...
                           vec_one,rhs);

% solve with Gauss Elimination
tol = 1e-6;

[Phi,iter, ~] = choose_solver(TotNumDM,2,tol,(NN+1)*TotNumDM,Phi_ini,...
                Fprolong,BlockL,DegDM,SupMat,SubMat);



end
