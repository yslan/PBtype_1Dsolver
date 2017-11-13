clear all; close all;
% This script is designed to run something outside the main.m
% Ex: loop of main with different parameter (differential capacitance)

% Initialize Test Case    % cases/
case_name   = 'loop_test'; %   <case_name>/
rea_name    = 'loop';  %      <rea_name>.rea
loop_name   = 'test';  %      <loop_name>.loop

method      = 2; % 1, 2, 3,...  mv to rea
% mode_sc     = 'fsolve'; % method to solve C in eq1 2
mode_sc     = 'ptws_NT_couple'; % method to solve C in eq1 2

% Include files and declare globals
addpath([pwd '/Poisson_1DMDMG/'])
init_loop_folder(case_name,rea_name,loop_name);
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
[loop,linfo] = read_loop(case_name,loop_name); % read loop


% Initials
[DegDM,Xprolong,JacVec,MassVec,Diff] = init_SEM(xmin,xmax,NN,TotNumDM);
%  initial condition
[C1,C2,Phi] = init_guess(Xprolong);


% Loop
for il = 1:linfo.Nloop
    % Display info
    disp(['Loop: il = ' num2str(il) '/' num2str(linfo.Nloop)])

    % Update global parameter(s)
    loop_uparam(loop,linfo,il);

    % Initial Guess
    %    the solution of "il-1"

    % Main solve
    [C1,C2,Phi] = main_solve(C1,C2,Phi,method,mode_sc);

    % savedata
    savedata(case_name,loop(il).name,il,Xprolong,C1,C2,Phi);
    % plotdata;
    % userchk;
end



%userchk;
%userplot;
%savedata; % save V0, Q, C_D
%plotdata;


% Remove path, allocation
rmpath([pwd '/Poisson_1DMDMG/'])
rmpath([pwd '/cases/' case_name])
