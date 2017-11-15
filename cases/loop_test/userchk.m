function [l2,linf] = userchk(varargin)
% check whatever you want in post process
% Input: sol_1,exact_1,sol_2,exact_2,...
global param

global TotNumDM NN xmin xmax
global k_B T e zval
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff

narginchk(2,8);
if mod(nargin,2)==1
    error('checking error must input in pairs');
end
Nchk = nargin/2;



f_sq = @(x)x.^2;


for i = 1:Nchk
    field_sol = varargin{2*i-1};
    field_exact = varargin{2*i};

    tmp = field_sol - field_exact;
    l2(i)   = (MassVec.*JacVec)' * f_sq(tmp);
    linf(i) = max(abs(tmp));
end

end
