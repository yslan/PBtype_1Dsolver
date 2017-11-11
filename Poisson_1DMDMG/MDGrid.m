% This script only creat the multidomain grid points

function [Xdomain, Jac] = MDGrid(DegDM,Xendpt);

% initialized required memory for variables
TotNumDM = length(Xendpt) - 1;

Nmax = max(DegDM); Nmaxp=Nmax+1;
    
Xdomain = zeros(Nmaxp,TotNumDM); Jac = Xdomain; 
    
% define the linear coordinate transformation funtion
T = @(x,k) Xendpt(1,k+1) + ((Xendpt(1,k+1) - Xendpt(1,k)) / 2) * (x - 1);
jaco = @(x,k) ((Xendpt(1,k+1) - Xendpt(1,k)) / 2);    

% loop all domains and compute X, Jac, Mass and D
for k = 1 : TotNumDM
    ND=DegDM(k); NDp=ND+1;  % get the deg N of domain k
    [x,~,~] = lglnodes(ND); % Compute LGL points and weights of Deg N
    
    % Compute Physical coordinates and Jacobian
    Xdomain(1:NDp,k) = T(-x,k); Jac(1:NDp,k) = jaco(-x,k);
    
end

