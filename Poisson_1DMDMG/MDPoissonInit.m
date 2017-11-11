function [BlockL,SupMat,SubMat,sigma_m,sigma_p] = ...
      MDPoissonInit(DegDM,Mass,D,Jac,A,AlphaBeta)

% initialized required memory for variables
TotNumDM = length(DegDM); Nmax = max(DegDM); Nmaxp=Nmax+1;
     
Ahat = zeros([Nmaxp Nmaxp TotNumDM]); 
BlockL=Ahat; SubMat=Ahat; SupMat=Ahat;

%Compute physical parameters
%[AlphaBeta,gl,gr,A,f] = PDEparas(TotNumDM,DegDM,Xdomain);

% Initialize em ep Im Ip
em=zeros(Nmaxp,TotNumDM); ep=zeros(Nmaxp,TotNumDM);
Im=zeros([Nmaxp Nmaxp TotNumDM]); Ip=Im;

for k = 1 : TotNumDM
    ND=DegDM(k); NDp=ND+1; 
    
    % initialize e_+ and e_- vectors
    em(1:NDp,k) = [1 zeros(1,ND)]'; 
    ep(1:NDp,k) = [zeros(1,ND) 1]';
      
    %Compute I_+ and I_- matrices
    Im(1:NDp,1:NDp,k) = em(1:NDp,k) * em(1:NDp,k)';
    Ip(1:NDp,1:NDp,k) = ep(1:NDp,k) * ep(1:NDp,k)';
    
end
 
%unpacked the boundary data 
alpha_m = AlphaBeta(1,1);  beta_m = AlphaBeta(1,2);
alpha_p = AlphaBeta(2,1);  beta_p = AlphaBeta(2,2);
     
% Compute the Ahat = A./J vector  
for k=1:TotNumDM
    NDp=DegDM(k)+1;  Ahat(1:NDp,1:NDp,k) = diag(A(1:NDp,k)./Jac(1:NDp,k));
end
    
% assign interface parameters
c_sigma = 4; % interface
c_tau = 4;   % left/right
betahat = zeros(2,TotNumDM); sigmahat = betahat; sigmabar = betahat;
sigma_m = zeros(DegDM(    1   )+1,1); 
sigma_p = zeros(DegDM(TotNumDM)+1,1); 
% fprintf('c_sigma = %f, c_tau = %f  function: MDPoissonInit\n',c_sigma, c_tau) %print

%boundary penalty
%left boundary
k = 1; ND = DegDM(k); NDp = ND+1;
if beta_m ~= 0 %Robin and Neumann boundary condition

   sigma_m(1:NDp) =  A(1,k) / beta_m / Mass(1,k) * em(1:NDp,k);                           

elseif beta_m == 0 %Dirichelt boundary condition

   tmp = Jac(1,k) .* Mass(1:NDp,k);         
   tau_m =  A(1,k) * diag(1./tmp) * D(1:NDp,1:NDp,k)' * em(1:NDp,k);
   tautil_m = c_tau * A(1,k) / Jac(1,k) / ( Mass(1,k)^2 ); 
   sigma_m(1:NDp) = tau_m(1:NDp) + tautil_m * em(1:NDp,k);
   
end

%right boundary
k = TotNumDM; ND = DegDM(k); NDp = ND+1;
if beta_p ~= 0  %Robin and Neumann boundary condition
   
   sigma_p(1:NDp) = A(NDp,k) / beta_p / Mass(NDp,k) * ep(1:NDp,k);
      
elseif beta_p == 0 %Dirichlet boundary condition

   tmp = Jac(NDp,k) .* Mass(1:NDp,k);      
   tau_p = A(NDp,k) * diag(1./tmp) * D(1:NDp,1:NDp,k)' * ep(1:NDp,k);
   tautil_p = c_tau * A(NDp,k) / Jac(NDp,k) / ( Mass( 1,k)^2 ) ; 
   sigma_p(1:NDp) = - tau_p(1:NDp) + tautil_p * ep(1:NDp,k); 

end

% assign parameters for imposing penalty BCs
% assign betahat at the left- and right-most boundaries
k = 1;                        betahat(1,k) = beta_m ./ Jac(1,k);
k = TotNumDM; NDp=DegDM(k)+1; betahat(2,k) = beta_p ./ Jac(NDp,k);

% assign parameters for internal BCs
if TotNumDM >1
    % assign interface boundary parameters
    k = 1; ND=DegDM(k); NDp = ND+1;
    betahat(2,k) = A(NDp,k) ./Jac(NDp,k);
    sigmahat(2,k) = 1 / ( 2 * Mass(1,k) );
    sigmabar(2,k) = c_sigma * betahat(2,k) / (4 * Mass(1,k).^2);
    
    for k = 2:TotNumDM-1
        ND = DegDM(k); NDp=ND+1;
        betahat(1:2,k)  = A(1:ND:NDp,k) ./ Jac(1:ND:NDp,k);
        sigmahat(1:2,k) = 1 / ( 2 * Mass(1,k) );
        sigmabar(1:2,k) = c_sigma * betahat(1:2,k) / (4 * Mass(1,k).^2);
    end
    
    k = TotNumDM; ND=DegDM(k); NDp = ND+1;
    betahat(1,k) = A(1,k) ./Jac(1,k);
    sigmahat(1,k) = 1 / ( 2 * Mass(1,k) );
    sigmabar(1,k) = c_sigma * betahat(1,k) / (4 * Mass(1,k).^2);
    
    %interface penalty
    k = 1; ND=DegDM(k); NDp=ND+1;
    
    sigmatil(1:NDp,2,k) =((-betahat(2,k) / 2) * diag(1./Mass(1:NDp,k))...
        * D(1:NDp,1:NDp,k)' * ep(1:NDp,k))...
        + (sigmabar(1,k+1) * Mass(1,k+1)/Mass(NDp,k) * ep(1:NDp,k));
    
    for k = 2:TotNumDM-1
        NDp = DegDM(k)+1;
        
        sigmatil(1:NDp,1,k) = ((betahat(1,k) / 2) * diag(1./Mass(1:NDp,k))...
            * D(1:NDp,1:NDp,k)' * em(1:NDp,k))...
            + (sigmabar(2,k-1) * Mass(1,k-1)/Mass( 1 ,k) * em(1:NDp,k));
        
        sigmatil(1:NDp,2,k) =((-betahat(2,k) / 2) * diag(1./Mass(1:NDp,k))...
            * D(1:NDp,1:NDp,k)' * ep(1:NDp,k))...
            + (sigmabar(1,k+1) * Mass(1,k+1)/Mass(NDp,k) * ep(1:NDp,k));
    end
    
    k=TotNumDM; NDp=DegDM(k)+1;
    sigmatil(1:NDp,1,k) = ((betahat(1,k) / 2) * diag(1./Mass(1:NDp,k))...
        * D(1:NDp,1:NDp,k)' * em(1:NDp,k))...
        + (sigmabar(2,k-1) * Mass(1,k-1)/Mass( 1 ,k) * em(1:NDp,k));
    
end % if TotNumDM >1
       
%--construct multidomian ML operator
%----block diagonal part of ML without boundary operators
for k = 1 : TotNumDM
    NDp = DegDM(k)+1; tmp = D(1:NDp,1:NDp,k);
    BlockL(1:NDp,1:NDp,k) = tmp * Ahat(1:NDp,1:NDp,k) * tmp;
end

%----diagonal part with left most boundary operators
k = 1; ND=DegDM(k); NDp=ND+1;
BlockL(1:NDp,1:NDp,k) = BlockL(1:NDp,1:NDp,k)...
       - sigma_m(1:NDp) * ( alpha_m * em(1:NDp,k)'...
       - betahat(1,k) * em(1:NDp,k)' * D(1:NDp,1:NDp,k) );   

%----diagonal part with right most boundary operators
k = TotNumDM; ND=DegDM(k); NDp=ND+1;
BlockL(1:NDp,1:NDp,k) = BlockL(1:NDp,1:NDp,k)...
       - sigma_p(1:NDp) * ( alpha_p * ep(1:NDp,k)'...
       + betahat(2,k) * ep(1:NDp,k)' * D(1:NDp,1:NDp,k) );
   
% --- block diagonal part with internal boundary operators
if TotNumDM >1
    %Add interface penalty on Block 1
    k=1; ND=DegDM(k); NDp = ND+1; Iploc = Ip(1:NDp,1:NDp,k);
    
    BlockL(1:NDp,1:NDp,1) = BlockL(1:NDp,1:NDp,1) ...
        - sigmatil(1:NDp,2,1) * ep(1:NDp,1)'...
        - sigmahat(2,1) * betahat(2,1) * Iploc * D(1:NDp,1:NDp,1)...
        - sigmabar(2,1) * Iploc;
    
    %Add interface penalty on inter Block
    for k = 2 : TotNumDM-1
        NDp = DegDM(k)+1; 
        Iploc=Ip(1:NDp,1:NDp,k); Imloc=Im(1:NDp,1:NDp,k);
        
        BlockL(1:NDp,1:NDp,k) = BlockL(1:NDp,1:NDp,k)...
            - sigmatil(1:NDp,2,k) * ep(1:NDp,k)'...
            - sigmahat(2,k) * betahat(2,k) * Iploc * D(1:NDp,1:NDp,k)...
            - sigmabar(2,k) * Iploc;
        
        BlockL(1:NDp,1:NDp,k) = BlockL(1:NDp,1:NDp,k)...
            - sigmatil(1:NDp,1,k) * em(1:NDp,k)'...
            + sigmahat(1,k) * betahat(1,k) * Imloc * D(1:NDp,1:NDp,k)...
            - sigmabar(1,k) * Imloc;
        
    end
    
    %Add interface penalty on Block end
    k= TotNumDM; ND=DegDM(k); NDp=ND+1;
    
    BlockL(1:NDp,1:NDp,k) = BlockL(1:NDp,1:NDp,k)...
        - sigmatil(1:NDp,1,k) * em(1:NDp,k)'...
        + sigmahat(1,k) * betahat(1,k) ...
        * em(1:NDp,k) * em(1:NDp,k)' * D(1:NDp,1:NDp,k)...
        - sigmabar(1,k) * Im(1:NDp,1:NDp,k);
    
    %Add submatrix and supmatrix
    for k = 1 : TotNumDM-1
        NDp =DegDM(k)+1; RNDp=DegDM(k+1)+1;
        
        SupMat(1:NDp,1:RNDp,k)...
            = sigmatil(1:NDp,2,k) * em(1:RNDp,k+1)'...
            + sigmahat(2,k) * betahat(1,k+1)...
            * ep(1:NDp,k) * em(1:RNDp,k+1)' * D(1:RNDp,1:RNDp,k+1)...
            + sigmabar(2,k) * ep(1:NDp,k) * em(1:RNDp,k+1)';
        
    end
end % If TotNumDM > 1

    
% Construct the tri-diagonal block array            
for k=1:TotNumDM-1
    NDp = DegDM(k) + 1;      RNDp = DegDM(k+1) + 1; 
    tmp1 = repmat(Mass(1:NDp,k),1, NDp);
    tmp2 = repmat(Mass(1:NDp,k),1,RNDp);
    BlockL(1: NDp,1: NDp,k) = - tmp1 .* BlockL(1:NDp,1:NDp,k);
    SupMat(1: NDp,1:RNDp,k) = - tmp2 .* SupMat(1:NDp,1:RNDp,k);
    SubMat(1:RNDp,1: NDp,k) = transpose(SupMat(1:NDp,1:RNDp,k));
    
end
    k = TotNumDM; ND=DegDM(k); NDp =ND+1;
    BlockL(1:NDp,1:NDp,k) = -diag(Mass(1:NDp,k))*BlockL(1:NDp,1:NDp,k);

%Obtain the numerical solution by Gauss Elimination
%[v]=GaussEliminationMD(TotNumDM,DegDM,F,BlockL,SubMat,SupMat);


end


