function [Xprolong,DegDM,BlockL,SupMat,SubMat,Fprolong] = ...
                get_L(xmin,xmax,NN,TotNumDM,u,dudx,a_fun,f_fun)
% user specifies boundary parameters
alpha_m = 1; beta_m = 0;  % specifies here
alpha_p = 1; beta_p = 0;  % specifies here

ApproachType = 'CG';

switch ApproachType
    case 'CG'
        Xendpt   = xmin + (xmax-xmin)/TotNumDM*(0:TotNumDM);
        
        for nn = 1: length(NN);
            
            
            DegDM=NN(nn) * ones([1 TotNumDM]);
            
            % Initialize Mass and Diff Mat
            [Mass, Diff] = MDMassDiff(DegDM);
            
            % Initialize Phyical grid points and transformation jacobian
            [Xdomain, Jac] = MDGrid(DegDM,Xendpt);
            
            % Initialize Poisson PDE parameters
            [AlphaBeta,gl,gr,A,f] = PoissonPara_MG(DegDM,Xdomain,...
                    alpha_m, beta_m, alpha_p, beta_p, u, dudx, a_fun, f_fun);
            
                % Construct ML operator in block-tri-diag form and
            % penalty parameters
            [BlockL, SupMat, SubMat, sigma_m, sigma_p] = ...
                MDPoissonInit(DegDM,Mass,Diff,Jac,A,AlphaBeta);
            
            % Initialize Gauss Solver
%             [BlockLInv,SupMat_GE]=MDGaussInit(DegDM,BlockL,SubMat,SupMat);
            
            % Compute right hand side vector
            F = RHSVec(DegDM,Mass,Jac,f,gl,gr,sigma_m,sigma_p);
            
            % Store RHS F vector in prolong vector form
            Fprolong = StackVec(DegDM,F);
            Xprolong = StackVec(DegDM,Xdomain);


            

            
        end

end

end
