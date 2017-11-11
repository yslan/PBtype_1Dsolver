function param2global
global param

global TotNumDM NN xmin xmax
global k_B T e zval
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

TotNumDM    = param(4);
NN          = param(5);

xmin        = param(6);
xmax        = param(7);
zval        = [xmin;xmax];


k_B         = param(8);
T           = param(9);
e           = param(10);

zval(1)     = param(11);
zval(2)     = param(12);


c0          = param(13);
phi_L       = param(14);


tol_pot     = param(15);
tol_c       = param(16);
tol_res     = param(17);

Psi_1 = 0; Psi_2 = 0;
res = -res_Concen_eq(0,[c0;c0]);
Psi_1 = res(1); Psi_2 = res(2);


end
