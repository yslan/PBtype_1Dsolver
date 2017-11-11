
clear all
close all
%global DegDM
%global BlockL
%global SupMat
%global SubMat

%epsilon = 1/10/2;
%u_ex = @(x) - (2*x-3).^(1/epsilon) + 1;
u_ex = @(x) 1+cos(pi*x);
xmin = -1; 
xmax = 1;


%in put
NN = (4:4:32);
errmax = zeros(size(NN));
ConvRate = zeros(size(NN));

%ApproachType = 'One';  
%ApproachType = 'Separate';
ApproachType = 'CG'

switch ApproachType
    case 'Separate'
        tic 
        TotNumDM = 8;
        Xendpt   = xmin + (xmax-xmin)/TotNumDM*(0:TotNumDM);
        
        for nn = 1: length(NN);
            
            
            DegDM=NN(nn) * ones([1 TotNumDM]);
            
            % Initialize Mass and Diff Mat
            [Mass, Diff] = MDMassDiff(DegDM);
            
            % Initialize Phyical grid points and transformation jacobian
            [Xdomain, Jac] = MDGrid(DegDM,Xendpt);
            
            % Initialize Poisson PDE parameters
            [AlphaBeta,gl,gr,A,f] = PoissonPara(DegDM,Xdomain);
            
            % Construct ML operator in block-tri-diag form and
            % penalty parameters
            [BlockL, SupMat, SubMat, sigma_m, sigma_p] = ...
                MDPoissonInit(DegDM,Mass,Diff,Jac,A,AlphaBeta);
            
            % Initialize Gauss Solver
            [BlockLInv,SupMat_GE]=MDGaussInit(DegDM,BlockL,SubMat,SupMat);
            
            % Compute right hand side vector
            F = RHSVec(DegDM,Mass,Jac,f,gl,gr,sigma_m,sigma_p);
            
            % solve
            v3 = MDGaussSol2(DegDM,F,BlockLInv,SupMat_GE,SubMat);
            
            Xprolong = StackVec(DegDM,Xdomain);
            Vprolong = StackVec(DegDM,v3);
            errmax(nn) = max(abs(Vprolong - u_ex(Xprolong)));
            
        end
        toc
        % compute convergence rate
        LNN = length(NN);
        ConvRate(2:LNN) = log(errmax(2:LNN)./errmax(1:LNN-1))...
            ./log(NN(1:LNN-1)./NN(2:LNN));
        
        for nn = 1:length(NN)
            fprintf(' %3d  %.4e  %.2f \n', NN(nn),errmax(nn), ConvRate(nn))
        end
        
    case 'CG'
        tic 
        TotNumDM = 4;
        Xendpt   = xmin + (xmax-xmin)/TotNumDM*(0:TotNumDM);
        
        for nn = 1: length(NN);
            
            
            DegDM=NN(nn) * ones([1 TotNumDM]);
            
            % Initialize Mass and Diff Mat
            [Mass, Diff] = MDMassDiff(DegDM);
            
            % Initialize Phyical grid points and transformation jacobian
            [Xdomain, Jac] = MDGrid(DegDM,Xendpt);
            
            % Initialize Poisson PDE parameters
            [AlphaBeta,gl,gr,A,f] = PoissonPara(DegDM,Xdomain);
            
            % Construct ML operator in block-tri-diag form and
            % penalty parameters
            [BlockL, SupMat, SubMat, sigma_m, sigma_p] = ...
                MDPoissonInit(DegDM,Mass,Diff,Jac,A,AlphaBeta);
            
            % Initialize Gauss Solver
            [BlockLInv,SupMat_GE]=MDGaussInit(DegDM,BlockL,SubMat,SupMat);
            
            % Compute right hand side vector
            F = RHSVec(DegDM,Mass,Jac,f,gl,gr,sigma_m,sigma_p);
            
            % Store RHS F vector in prolong vector form
            Fprolong = StackVec(DegDM,F);
            F2 = vec2mdarray(DegDM,Fprolong,Xdomain);
            
            % solve by pcg function 
            %n2 = 21;
            %b2 = applyMoler(ones(n2,1));
            tol = 1e-15;
            maxit = 2000;
            %M2 = spdiags((1:n2)',0,n2,n2);
            %vlong = ones(size(Fprolong));
            
%            [vlong,flag,relres,iter,resvec] = ...
%                    pcg(@(x)MLx(x,DegDM,BlockL,SupMat,SubMat),Fprolong,...
%                        tol,maxit);  % by matlab pcg
            x0 = zeros(size(Fprolong)); 
            % x0 = StackVec(DegDM,Xdomain); %zeros(size(Fprolong));
            %xlong = StackVec(DegDM,Xdomain);
            %x0 = u_ex(xlong) + 1E-10*sin(xlong);
            %C = Mass; 
            C = ones(size(Mass));
            [vlong, iter1, resvec1,viter1] = ...
                cgp(x0, @(x)MLx(x,DegDM,BlockL,SupMat,SubMat), ...
                @(x)bbC(x,DegDM,C), Fprolong, maxit, tol);
            
            x0 = zeros(size(Fprolong));
            %x0 = StackVec(DegDM,Xdomain); %
            %xlong = StackVec(DegDM,Xdomain);
            %x0 = (2*xlong-3)+1; %u_ex(xlong) + 1E-6*sin(xlong);
            for k =1:TotNumDM
                NDp=DegDM(k)+1; 
                C(1:NDp,k) = 1./Mass(1:NDp,k);
            end
            %C = Mass.^(+1);
            %C = ones(size(Mass));
            [v2long, iter2, resvec2,viter2] = ...
                cgp(x0, @(x)MLx(x,DegDM,BlockL,SupMat,SubMat), ...
                @(x)bbC(x,DegDM,C), Fprolong, maxit, tol);        
            
            v2 = vec2mdarray(DegDM,v2long,Xdomain);

            v3 = MDGaussSol2(DegDM,F,BlockLInv,SupMat_GE,SubMat);
            
            Xprolong = StackVec(DegDM,Xdomain);
            Vprolong = StackVec(DegDM,v2);
            errmax(nn) = max(abs(Vprolong - u_ex(Xprolong)));
            
        end
        toc
        % compute convergence rate
        LNN = length(NN);
        ConvRate(2:LNN) = log(errmax(2:LNN)./errmax(1:LNN-1))...
            ./log(NN(1:LNN-1)./NN(2:LNN));
        
        for nn = 1:length(NN)
            fprintf(' %3d  %.4e  %.2f \n', NN(nn),errmax(nn), ConvRate(nn))
        end

    case 'One'
        tic
        TotNumDM = 8;
        Xendpt   = xmin + (xmax-xmin)/TotNumDM*(0:TotNumDM);
        
        for nn = 1: length(NN);
            DegDM=NN(nn) * ones([1 TotNumDM]);
          
            % Solve Poisson problem
            [Xdomain,v] = MDPoissonSolver(Xendpt,DegDM);
            Xprolong = StackVec(DegDM,Xdomain);
            Vprolong = StackVec(DegDM,v);
            errmax(nn) = max(abs(Vprolong - u_ex(Xprolong)));
        end
        toc
        % compute convergence rate
        LNN = length(NN);
        ConvRate(2:LNN) = log(errmax(2:LNN)./errmax(1:LNN-1))...
            ./log(NN(1:LNN-1)./NN(2:LNN));
        
        for nn = 1:length(NN)
            fprintf(' %3d  %.4e  %.2f \n', NN(nn),errmax(nn), ConvRate(nn))
        end
        
end
figure(1)
plot(Xprolong,u_ex(Xprolong),Xprolong,Vprolong,'.','MarkerSize',14)
axis([xmin xmax min(Vprolong) max(Vprolong)])

figure(2)
loglog((0:iter1),resvec1(1:iter1+1),'o-',...
       (0:iter2),resvec2(1:iter2+1),'*-')
   
for kk=1:iter1
    erriter1(kk) = norm(abs(u_ex(Xprolong) - viter1(:,kk)));
end 
for kk=1:iter2
    erriter2(kk) = norm(abs(u_ex(Xprolong) - viter2(:,kk)));
end 
figure(3)
loglog((1:iter1),erriter1(1:iter1),'o-',...
       (1:iter2),erriter2(1:iter2),'*-')
   axis([1 3000 1E-16 1E2])
   
figure(4)
for kk=1:1:iter2
    plot(Xprolong,viter2(:,kk),'o')
    hold on
end

% for nn = 1:length(V)
%     fprintf(' %.15e  %.15e \n', Xprolong(nn),V(nn))
% end
