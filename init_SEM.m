function [DegDM,Xprolong,Jac_vec,Mass_vec,Diff] = init_SEM(xmin,xmax,NN,TotNumDM)

Xendpt   = xmin + (xmax-xmin)/TotNumDM*(0:TotNumDM);

for nn = 1: length(NN);
	DegDM=NN(nn) * ones([1 TotNumDM]);
            
	% Initialize Mass and Diff Mat
	[Mass, Diff] = MDMassDiff(DegDM);
	Mass_vec = reshape(Mass,TotNumDM*(NN+1),1);
            
    % Initialize Phyical grid points and transformation jacobian
    [Xdomain, Jac] = MDGrid(DegDM,Xendpt);
    Xprolong = StackVec(DegDM,Xdomain);
    Jac_vec = reshape(Jac,(NN+1)*TotNumDM,1);
end