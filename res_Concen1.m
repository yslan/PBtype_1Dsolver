function res=res_Concen1(phi,C1,C2_ini)
global k_B T e zval
global Psi_1 Psi_2
[f1,f2] = userf(C1,C2_ini);
res = Psi_1 - k_B * T * log(C1) - e*zval(1)*phi - f1;

end