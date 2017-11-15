function userchk(C1,C2,Phi)
% check whatever you want in post process
global param

global TotNumDM NN xmin xmax
global k_B T e zval
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff

Q    = compute_Q(C1,C2,0);
Q_ex = compute_Q(C1,C2,1);


Q_err = abs(Q - Q_ex);

disp(['    Q = ' num2str(Q) ', Q_ex = ' num2str(Q_ex) ', Q_err = ' num2str(Q_err)])
% f_sq = @(x)x.^2;

% Nchk = 1;
% for i = 1:Nchk
%     field_sol = varargin{2*i-1};
%     field_exact = varargin{2*i};
% 
%     tmp = field_sol - field_exact;
%     l2(i)   = (MassVec.*JacVec)' * f_sq(tmp);
%     linf(i) = max(abs(tmp));
% end
    
end
