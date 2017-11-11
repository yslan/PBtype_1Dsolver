function res=res_Concen_eq(phi,C)
% Given C1, C1 and Phi, compute res in eq Concent 1 and Concent 2

global k_B T e zval
global Psi_1 Psi_2
C = abs(C);

Nx = 1/2*length(C);
C1 = C(1:Nx);
C2 = C(Nx+1:2*Nx);

[f1,f2] = userf(C1,C2);
res1 = Psi_1 - k_B * T * log(C1) - e*zval(1)*phi - f1;
res2 = Psi_2 - k_B * T * log(C2) - e*zval(2)*phi - f2;
res = [res1;res2];
end
