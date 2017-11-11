function [AlphaBeta, g_m, g_p, Adomain, fdomain] = ...
    PoissonPara_MG(DegDM,Xdomain, ...
    alpha_m, beta_m, alpha_p, beta_p, ...
    u, dudx, a_fun, f_fun)
% -------------------------------------------------------------------------     
% This script is used to specify parameters of the general Poission problem
%
% PDE : -(a(x)*u'(x))' = f(x), x \in [x_left, x_right]
% BCs : alpha_m * u(x_left)  - beta_m * u'(x_left ) = g_m
%       alpha_p * u(x_right) + beta_p * u'(x_right) = g_p
% 
% u(x): unknown solution
% a(x): piecewise continuos function
% 
% function inputs ---------------------------------------------------------
%  TotNumDM:      size (1,1)
%  ^^^^^^^^       Total number of domains
% 
%  DegDM:         size (1,TotNumDM) 
%  ^^^^^          Degree of the approximation polynomials in each domain
% 
%  Xdomain:       size (Nmax+1,TotNumDM)
%  ^^^^^^^        coordinates of the mapped LGL grid points in each domain
% -------------------------------------------------------------------------
% function outputs
% AlphaBeta:       size(2,2). Need to specifies in this scripts (see below)
% ^^^^^^^^^        parameters for setting boundary conditions.
%                  AlphaBeta = [alpha_m beta_m; alpha_p beta_m]
%
% g_m & g_p:       size (1x1) values of the boundary conditions
% ^^^^^^^^^
%
% Adomain & f:     size(Nmax+1,TotNumDM)
% ^^^^^^^^^^^      pointswise values of a(x) and f(x) at the grid points 
%                  stored based on domains
% 
% Aplong & fplong: size(Nmax+1,TotNumDM)
% ^^^^^^^^^^^^^^^  pointswise values of a(x) and f(x) at the grid points
%                  stacked as prolong vectors
%
%-------------------------------------------------------------------------

TotNumDM = length(DegDM);





% function begin
AlphaBeta = [alpha_m beta_m; alpha_p beta_p];


% compute the pointwise values of a(x)
Adomain = zeros(size(Xdomain)); 
for DDK = 1:TotNumDM
    ND = DegDM(DDK); NDp = ND+1;
    Adomain(1:NDp,DDK) = a_fun(Xdomain(1:NDp,DDK)); 
end

% compute the pointwise values of f(x)
fdomain= zeros(size(Xdomain));
for DDK = 1:TotNumDM
    ND = DegDM(DDK); NDp = ND+1;
    fdomain(1:NDp,DDK) = f_fun(Xdomain(1:NDp,DDK));
end

% assign g_l and g_r 
xm = Xdomain(1,1);   % domain left-end point 
DDK = TotNumDM;  NDp = DegDM(DDK)+1;
xp = Xdomain(NDp,DDK); % domain right-end point

g_m = alpha_m * u(xm) - beta_m *dudx(xm);
g_p = alpha_p * u(xp) + beta_p *dudx(xp);


end

