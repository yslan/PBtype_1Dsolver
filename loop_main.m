clear all; close all;clc;format shorte;
warning('off');rmpath(genpath([pwd '/cases']));warning('on');
% This script is designed to run something outside the main.m
% Ex: loop of main with different parameter (differential capacitance)

% Initialize Test Case    % cases/
case_name   = 'NaF_PBbikerman_loop'; %   <case_name>/
rea_name    = 'naf';  %      <rea_name>.rea
% loop_name   = 'test_mode2';  %      <loop_name>.loop
loop_name   = '0p5';  %      <loop_name>.loop

method      = 2; % 1, 2, 3,...  mv to rea
% mode_sc     = 'fsolve'; % method to solve C in eq1 2
mode_sc     = 'ptws_NT_couple'; % method to solve C in eq1 2

% Include files and declare globals
addpath([pwd '/Poisson_1DMDMG/'])
init_loop_folder(case_name,rea_name,loop_name);
addpath([pwd '/cases/' case_name]) % incluse param.m and usr.m

global param

global TotNumDM NN xmin xmax
global k_B T e zval epsilon
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff

global ifsol

% Input parameters
tic
read_rea(case_name,rea_name); % read param
param2global; % set global param
[loop,linfo] = read_loop(case_name,loop_name); % read loop
t0 = toc;
disp('Read input...done!!')

% Initials
[DegDM,Xprolong,JacVec,MassVec,Diff] = init_SEM(xmin,xmax,NN,TotNumDM);
%  initial condition
[C1,C2,Phi] = init_guess(Xprolong);
t1 = toc;
disp('Initialization ...done!!')


tbl = toc;
% Loop
for il = 1:linfo.Nloop
    tl0 = toc;
    % Display info
    disp(['Loop: il = ' num2str(il) '/' num2str(linfo.Nloop)])

    % Update global parameter(s)
    loop_uparam(loop,linfo,il);

    % Initial Guess
    %    the solution of "il-1"


    % Main solve
    [C1,C2,Phi] = main_solve(C1,C2,Phi,method,mode_sc);
    tl1 = toc;
    disp(['    Loop: il = ' num2str(il) ' is solved, cost ' num2str(tl1-tl0) 'sec'])


    % savedata

    savedata(case_name,loop(il).name,il,Xprolong,C1,C2,Phi);
    if (ifsol); Q(il) = userchk(C1,C2,Phi); end;
    disp('    date saved.')

    % plotdata;
%    [l2,linf] = userchk(C1,C2,Phi);
    tl2 = toc;
end
tal = toc;


if (ifsol); [V0,C_D] = userfinalchk(Q,loop,linfo);end

%userchk;
%userplot;
%savedata; % save V0, Q, C_D
%plotdata;

tfinal = toc;
sp10 = '          ';
disp('=================================================================================')
disp(['Case:         ' case_name])
disp(['rea:          ' rea_name])
disp(['  info:'])
disp(['    TotNumDM = ' num2str(TotNumDM) ',  NN = ' num2str(NN)])
disp(['    xmin     = ' num2str(xmin) ', xmax = ' num2str(xmax)])
disp(['    c0       = ' num2str(c0)])
disp(['    zval     = ' num2str(zval')])
disp(['Loop:         ' loop_name])
disp(['  info:'])
disp(linfo)
disp(['  Time: (total/average): ' sp10 sp10 num2str(tal-tbl) ' / ' num2str((tal-tbl)/linfo.Nloop)])
disp('=================================================================================')
disp(['                                                            Total time:' num2str(tfinal)])
disp('Delete dependnecy, program should exit any time...')
 
% Remove path, allocation
rmpath([pwd '/Poisson_1DMDMG/'])
rmpath([pwd '/cases/' case_name])
