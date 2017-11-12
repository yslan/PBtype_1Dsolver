clear all; close all;
% Solving PB-type Equation, which is steady state PNP type with 
% no-flux boundary condition
%  Eqs:
%   Psi_1 = k_B*T * ln(C_1) + e*z_1*phi + f_1(C_1,C_2) ... (1)
%   Psi_2 = k_B*T * ln(C_2) + e*z_2*phi + f_2(C_1,C_2) ... (2)
%   - (epsilon * Phi')' = e*(z_1*C_1+z_2*C_2) + e*Q    ... (3)


% Initialize Test Case
case_name   = 'PB_steric_g_1112';
rea_name    = 'test';

method      = 2; % 1, 2, 3,...  mv to rea
% mode_sc     = 'fsolve'; % method to solve C in eq1 2
mode_sc     = 'ptws_NT_couple'; % method to solve C in eq1 2

% Include files and declare
addpath([pwd '/Poisson_1DMDMG/'])
init_folder(case_name,rea_name);
addpath([pwd '/cases/' case_name]) % incluse param.m and usr.m

global param

global TotNumDM NN xmin xmax
global k_B T e zval
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong 
global MassVec JacVec Diff


% Input parameters
read_rea(case_name,rea_name); % read param
param2global; % set global param


% Initials
[DegDM,Xprolong,JacVec,MassVec,Diff] = init_SEM(xmin,xmax,NN,TotNumDM);

%  initial condition
Phi_ini = zeros(size(Xprolong));
C1_ini = ones(size(Xprolong));
C2_ini = ones(size(Xprolong));
C_ini = [C1_ini;C2_ini];
Phi = Phi_ini;
C1 = C1_ini;
C2 = C2_ini;


% Main solving part
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
        [C,fval] =fsolve(@(C)TotalRes_Concen(C,Phi_ini),C_ini);

        Nx = length(C)/2;
        C1 = C(1:Nx); C2 = C(Nx+1:2*Nx);
        Phi = solve_phi_by_c(C1,C2,Phi_ini);
end


figure(2)

plpo = plot(Xprolong,Phi,'k');hold on
plcN = plot(Xprolong,C1,'r');
plcP = plot(Xprolong,C2,'b');
legend([plcN plcP plpo],'cN','cP','Potent')
title('EDL simulation, single run, PB eq')


% Remove path, allocation
rmpath([pwd '/Poisson_1DMDMG/'])
rmpath([pwd '/cases/' case_name])

