function loop_uparam(loop,linfo,il)
% update parameters at each cases of loop

global param

global TotNumDM NN xmin xmax
global k_B T e zval
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff


lmode = linfo.mode;
Nvar  = linfo.Nvar;

switch lmode
    case 0
        for iv = 1:Nvar
            ip = loop(il).param_i(iv);
            param(ip) = loop(il).var_value(iv);
        end
        param2global;

    case 1
        for iv = 1:Nvar
            eval([loop(il).var_name{iv} '=' num2str(loop(il).var_value) ';']);
        end
        Psi_1 = 0; Psi_2 = 0;
        res = -res_Concen_eq(0,[c0;c0]);
        Psi_1 = res(1); Psi_2 = res(2);
    case 2 
        for iv = 1:Nvar
            eval(['varf = ' loop(il).var_eval{iv} ';'])
            eval([loop(il).var_name{iv} '=' num2str(varf(il)) ';']);
        end
        Psi_1 = 0; Psi_2 = 0;
        res = -res_Concen_eq(0,[c0;c0]);
        Psi_1 = res(1); Psi_2 = res(2);
end


end
