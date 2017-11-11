function [C1,C2,fval] = solve_c_by_phi(Phi,C1_ini,C2_ini,mode)
% Given Phi, get C1 and C2
%   mode: 
%       fsolve
%       ptws_NT: point wise Newton

        


Nx = length(Phi);
switch mode
    case 'fsolve'
        option = optimoptions('fsolve','Display','none');
        
        C_ini = [C1_ini;C2_ini];
        [C_out,fval] = fsolve(@(C_in)res_Concen_eq(Phi,C_in),C_ini,option);
        
        C1 = C_out(1:Nx); C2 = C_out(Nx+1:2*Nx);
    
    case 'ptws_NT'
    
end
