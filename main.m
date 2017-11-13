clear all; close all;
% Solving PB-type Equation, which is steady state PNP type with 
% no-flux boundary condition
%  Eqs:
%   Psi_1 = k_B*T * ln(C_1) + e*z_1*phi + f_1(C_1,C_2) ... (1)
%   Psi_2 = k_B*T * ln(C_2) + e*z_2*phi + f_2(C_1,C_2) ... (2)
%   - (epsilon * Phi')' = e*(z_1*C_1+z_2*C_2) + e*Q    ... (3)


% Initialize Test Case
case_name   = 'test1113';
rea_name    = 'rea1113';

method      = 2; % 1, 2, 3,...  mv to rea
% mode_sc     = 'fsolve'; % method to solve C in eq1 2
mode_sc     = 'ptws_NT_couple'; % method to solve C in eq1 2

% Include files and declare globals
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
%  initial guess
[C1,C2,Phi] = init_guess(Xprolong);


% Main solving part
[C1,C2,Phi] = main_solve(C1,C2,Phi,method,mode_sc);


% Save and plot
savedata(case_name,rea_name,0,Xprolong,C1,C2,Phi); % save into .csv
% plotdata;

figure(2)

plpo = plot(Xprolong,Phi,'k');hold on
plcN = plot(Xprolong,C1,'r');
plcP = plot(Xprolong,C2,'b');
legend([plcN plcP plpo],'cN','cP','Potent')
title('EDL simulation, single run, PB eq')


% Remove path, allocation
rmpath([pwd '/Poisson_1DMDMG/'])
rmpath([pwd '/cases/' case_name])

