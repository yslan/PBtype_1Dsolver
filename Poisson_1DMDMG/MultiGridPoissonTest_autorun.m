function MultiGridPoissonTest_autorun(in_filename,name_label)
% in_filename = '0219_example_1.txt';
fig_save = 1;
mat_save = 1;

rel_dir_in = 'input/';
rel_dir_out = 'output/';
t_flag = char(datetime('now','Format','yyyy_MMdd'));
mkdir(strcat(rel_dir_out,t_flag));
mkdir(strcat(rel_dir_out,t_flag,'/',name_label));

in_file_path = strcat(rel_dir_in,in_filename);
out_file_path = strcat(rel_dir_out,t_flag,'/',name_label,'/');

copyfile(in_file_path,[out_file_path in_filename])
% open files
fclose('all');

fid = fopen(in_file_path,'r');
if fid ==-1
    error('ERROR: cannot open the input file')
end
tline = fgets(fid); % no info. line
tline = fgets(fid); % no info. line
tline = fgets(fid); 
Test_num = sscanf(tline,'%d');

for i = 1:Test_num
    fprintf(['Test case: ' num2str(i) ' start...\n'])
    tline = fgets(fid); % no info. line
    tline = fgets(fid); % no info. line
    tline = fgets(fid); % no info. line
    tline = fgets(fid);
    vec = sscanf(tline,' %d %d %d');
    NN = vec(1); TotNumDM = vec(2); Nc = vec(3);
    
    tline = fgets(fid); % no info. line
    tline = fgets(fid);
    vec = sscanf(tline,' %f %f %d');
    xmin = vec(1); xmax = vec(2); fun_case = vec(3);
    
    tline = fgets(fid); % no info. line
    tline = fgets(fid);
    vec = sscanf(tline,' %d %d %d');
    sType_vc = vec(1); sType_fin = vec(2); algo_vc = vec(3);
    
    tline = fgets(fid); % no info. line
    tline = fgets(fid);
    vec = sscanf(tline,' %d %d %f');
    N_vcycle = vec(1); m_smooth = vec(2); sigma = vec(3);
    
    tline = fgets(fid); % no info. line
    tline = fgets(fid);
    vec = sscanf(tline,' %f %d');
    tol = vec(1); maxit = vec(2);
    
    % data transform
    save_label = [name_label '_test_' num2str(i)];
    
    save_info.fig_save = fig_save;
    save_info.mat_save = mat_save;
    save_info.out_file_path = out_file_path;
    save_info.save_label = save_label;
    
    parameter.NN = NN;
    parameter.TotNumDM = TotNumDM;
    parameter.Nc = Nc;
    parameter.xmin = xmin;
    parameter.xmax = xmax;
    parameter.fun_case = fun_case;
    parameter.sType_vc = sType_vc;
    parameter.sType_fin = sType_fin;
    parameter.algo_vc = algo_vc;
    parameter.tol = tol;
    parameter.maxit = maxit;
    parameter.N_vcycle = N_vcycle;
    parameter.m_smooth = m_smooth;
    parameter.sigma = sigma;

    MultiGridCG_function(parameter,save_info);
    fprintf(['Test case: ' num2str(i) '       done!!\n'])
end
% Parematers Explaination:
% Info of Domain
%     NN = 20;
%     TotNumDM = 256;
%     Nc = 4;%floor(NN/2); % coarse grid size
%     xmin = 1; 
%     xmax = 2.0;
% Choose linear solver, solving Ax=b
%   inside the v-cycle
%     sType_vc = 1;
%   after the v-cycle
%     sType_fin = 1;
%   Table of indices of solvers
%     Ind   Solver              Abbr. in code       Domain type
%       0   MATLAB slash        mat_slash           single
%       1   CG                  CG                  single/Multi
%       2   MATLAB gmres        mat_gmres           single/Multi
%       3   MATLAB pcg          mat_pcg             single/Multi
%       4   GaussElimination    Gauss               single/Multi
% needed in any Iterative method
%     tol = 1e-6;
%     maxit = 6000;
% Choose Algorithm in the smoothing part of V-clcle
%     algo_vc = 1; % (recommand choose 1)
%   Table of indices of Algorithm
%     Ind        Algorithm
%       0       Precondition Jacobi
%       1       Non-overlapped Jacobi Schwartz (actually overlapped by penalty)
%       2       Overlapped Jacobi Schwartz
% fun_case = 2; % chosing test case







% Parematers Explaination:
% Info of Domain
%     NN = 20;
%     TotNumDM = 256;
%     Nc = 4;%floor(NN/2); % coarse grid size
%     xmin = 1; 
%     xmax = 2.0;
% Choose linear solver, solving Ax=b
%   inside the v-cycle
%     sType_vc = 1;
%   after the v-cycle
%     sType_fin = 1;
%   Table of indices of solvers
%     Ind   Solver              Abbr. in code       Domain type
%       0   MATLAB slash        mat_slash           single
%       1   CG                  CG                  single/Multi
%       2   MATLAB gmres        mat_gmres           single/Multi
%       3   MATLAB pcg          mat_pcg             single/Multi
%       4   GaussElimination    Gauss               single/Multi
% needed in any Iterative method
%     tol = 1e-6;
%     maxit = 6000;
% Choose Algorithm in the smoothing part of V-clcle
%     algo_vc = 1; % (recommand choose 1)
%   Table of indices of Algorithm
%     Ind        Algorithm
%       0       Precondition Jacobi
%       1       Non-overlapped Jacobi Schwartz (actually overlapped by penalty)
%       2       Overlapped Jacobi Schwartz
% fun_case = 2; % chosing test case


end