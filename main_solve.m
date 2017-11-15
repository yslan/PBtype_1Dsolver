function [C1,C2,Phi] = main_solve(C1_ini,C2_ini,Phi_ini,varargin)
% main solving part

global param

global TotNumDM NN xmin xmax
global k_B T e zval epsilon
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff

narginchk(4,5);
method = varargin{1};
if nargin == 5
   mode_sc = varargin{2};
end

switch method
    case 1 % relaxation
        alpha = 0.01;

        for i = 1:1000
            [C1,C2,fval] = solve_c_by_phi(Phi,C1,C2,mode_sc);

            Phi_new = solve_phi_by_c(C1,C2,Phi);
            Phi = Phi_new*alpha + (1-alpha)*Phi;

            figure(1)
            % plot(Xprolong,Phi,'k');hold on
            plot(Xprolong,C1);hold on
            plot(Xprolong,C2);hold off
            title([num2str(i) 'steps, a.m. fval = ' num2str(max(abs(fval)))])
            drawnow
        end
    case 2 % fsolve on Phi
        [Phi,fval] =fsolve(@(phi)TotalRes_Poisson(phi,C1_ini,C2_ini),Phi_ini);
        [C1,C2,fval] = solve_c_by_phi(Phi,C1_ini,C2_ini,mode_sc);

    case 3 % fsolve on C
        [C,fval] =fsolve(@(C)TotalRes_Concen(C,Phi_ini),[C1_ini;C2_ini]);

        Nx = length(C)/2;
        C1 = C(1:Nx); C2 = C(Nx+1:2*Nx);
        Phi = solve_phi_by_c(C1,C2,Phi_ini);
end


end
