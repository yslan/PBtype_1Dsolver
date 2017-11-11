function res=res_Concen2(phi,C2,C1_ini)
global k_B T e zval
global Psi_1 Psi_2
[f1,f2] = userf(C1_ini,C2);

res = Psi_2 - k_B * T * log(C2) - e*zval(2)*phi - f2;
end