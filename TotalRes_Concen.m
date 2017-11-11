function res=TotalRes_Concen(C,Phi_ini)
% Given C1, C2, compute total res of PDEin form of Concen 1 and 2

global k_B T e zval
global Psi_1 Psi_2
C = abs(C);

Nx = 1/2*length(C);
C1 = C(1:Nx);
C2 = C(Nx+1:2*Nx);

Phi = solve_phi_by_c(C1,C2,Phi_ini);

[f1,f2] = userf(C1,C2);
res1 = Psi_1 - k_B * T * log(C1) - e*zval(1)*Phi - f1;
res2 = Psi_2 - k_B * T * log(C2) - e*zval(2)*Phi - f2;
res = [res1;res2];
end
