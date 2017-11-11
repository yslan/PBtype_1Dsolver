function [C1,fval] = solve_c1_by_phi(Phi,C1_ini,C2_ini)

% C1 = zeros(Nx,1);
% C2 = zeros(Nx,1);
option = optimoptions('fsolve','Display','none');
% c_ini = [C1_ini;C2_ini];

% for i = 1:length(Phi)
%     [c_out,fval] = fsolve(@(c_in)res_Concen(Phi(i),c_in(1),c_in(2)),c_ini,option);
% 
%     C1(i) = c_out(1); C2(i) = c_out(2);
% end
[C1,fval] = fsolve(@(c_in)res_Concen1(Phi,c_in,C2_ini),C1_ini,option);
end