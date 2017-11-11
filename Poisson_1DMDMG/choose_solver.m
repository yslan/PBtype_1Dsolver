function [x_sol, iter, label] = choose_solver(TotNumDM,ind,tol,maxit,x_ini,b,varargin)
                       % choose_solver(TotNumDM,ind,tol,maxit,x_vc_ini,b,A)
                       % choose_solver(TotNumDM,ind,tol,maxit,x_vc_ini,b,BlockL,DegDM,SupMat,SubMat)
% ======================================================================= %
%   Solving Ax=b, 
%   Table of indices of solvers
%     Ind   Solver              Abbr. in code       Domain type
%       0   MATLAB slash        mat_slash           single
%       1   CG                  CG                  single/Multi
%       2   MATLAB gmres        mat_gmres           single/Multi
%       3   MATLAB pcg          mat_pcg             single/Multi
%       4   GaussElimination    Gauss               single/Multi                   
% ======================================================================= %
narginchk(7, 10);

if nargin == 7
    A = varargin{1};
elseif nargin == 10
    BlockL = varargin{1};
    DegDM = varargin{2};
    SupMat = varargin{3};
    SubMat = varargin{4};
end

switch ind
    case 0   % MATLAB slash        mat_slash           single
        label = 'mat slash';
        iter = 0;
        x_sol = A\b;
        
    case 1   % CG                  CG                  single/Multi
        label = 'CG';
        
        % C = Mass; % pre-condition
        C = ones(length(x_ini),TotNumDM); % No pre-condition
        [x_sol, iter, ~] = ...
            cgp(x_ini, @(x)MLx(x,DegDM,BlockL,SupMat,SubMat), ...
            @(x)bbC(x,DegDM,C), b, maxit, tol);
        
    case 2   % MATLAB gmres        mat_gmres           single
        label = 'mat gmres';
        
        [x_sol,~,~,iter_gmres] =...
            gmres(@(x)MLx(x,DegDM,BlockL,SupMat,SubMat),b,[],tol,maxit,[],[],x_ini);
        iter = iter_gmres(2);
        
    case 3   % MATLAB pcg          mat_pcg             single
        label = 'mat pcg';
        [x_sol,~,~,iter] =...
            pcg(@(x)MLx(x,DegDM,BlockL,SupMat,SubMat),b,tol,maxit,[],[],x_ini);
        
    case 4   % GaussElimination    Gauss               single/Multi 
        label = 'Gauss';
        iter = 0;
        b_array = vec2mdarray(DegDM,b,BlockL);
        
        x_array = GaussEliminationMD(TotNumDM,DegDM,b_array,BlockL,SubMat,SupMat);
        x_sol = mdarray2vec(DegDM,x_array);

end










end