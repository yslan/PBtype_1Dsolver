function [Mass Diff] = MDMassDiff(DegDM)

% initialized required memory for variables
TotNumDM = length(DegDM); Nmax = max(DegDM); Nmaxp=Nmax+1;
    
Mass = zeros(Nmaxp,TotNumDM);  
Diff = zeros([Nmaxp Nmaxp TotNumDM]); 
    
% loop all domains and compute X, Jac, Mass and D
for k = 1 : TotNumDM
    ND=DegDM(k); NDp=ND+1;  % get the deg N of domain k
    [x,w,~] = lglnodes(ND); % Compute LGL points and weights of Deg N
    
    %Compute the mass differential matrix of Deg N
    Mass(1:NDp,k) = w;  Diff(1:NDp,1:NDp,k) = collocD(-x); 
end