function f = RHSVec(DegDM,Mass,Jac,f,gm,gp,sigma_m,sigma_p)

TotNumDM = length(DegDM); %F = zeros(size(f));
for k = 1:TotNumDM
    NDp = DegDM(k)+1;
    f(1:NDp,k) =  Mass(1:NDp,k) .* Jac(1:NDp,k) .* f(1:NDp,k);            
end
% Compute right hans side
    k = 1; NDp = DegDM(k)+1;  
    f(1:NDp,k) = f(1:NDp,k) + Mass(1:NDp,k) .* sigma_m(1:NDp) * gm;
    
    k = TotNumDM; NDp = DegDM(k)+1; 
    f(1:NDp,k) = f(1:NDp,k) + Mass(1:NDp,k) .* sigma_p(1:NDp) * gp;
end