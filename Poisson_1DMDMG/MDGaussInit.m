function  [BlockLInv,SupMat_GE]=MDGaussInit(DegDM,BlockL,SubMat,SupMat)
%Gauss Elimination 
%Forward sweeping
TotNumDM = length(DegDM); 
SupMat_GE = zeros(size(SupMat));
BlockLInv = zeros(size(BlockL));

switch TotNumDM
    case 1
        k = 1; ND = DegDM(k); NDp = ND+1;
        tmp = inv( BlockL(1:NDp,1:NDp,k)) ;
        BlockLInv(1:NDp,1:NDp,k) = tmp;
%        v(1:NDp,k) = BlockLInv(1:NDp,1:NDp) * F(1:NDp,k) ;
        
    otherwise
        
        k = 1; NDp=DegDM(k)+1; RNDp=DegDM(k+1)+1;
        tmp = inv( BlockL(1:NDp,1:NDp,k)) ;
        BlockLInv(1:NDp,1:NDp,k) = tmp;
        SupMat_GE(1:NDp,1:RNDp,k) = tmp * SupMat(1:NDp,1:RNDp,k) ;
%        F_GE(1:NDp,k) = tmp * F(1:NDp,k) ; % tmp = BlockLInv
        
        for k = 2 : TotNumDM-1
            LNDp = DegDM(k-1)+1; NDp = DegDM( k )+1; RNDp = DegDM(k+1)+1;
            
            tmp = inv( BlockL(1:NDp,1:NDp,k) ...
                - SubMat(1:NDp,1:LNDp,k-1) * SupMat_GE(1:LNDp,1:NDp,k-1) );
            BlockLInv(1:NDp,1:NDp,k) = tmp;
            SupMat_GE(1:NDp,1:RNDp,k) = tmp * SupMat(1:NDp,1:RNDp,k) ;
            
%            F_GE(1:NDp,k) = tmp * ( F(1:NDp,k) ...
%                - SubMat(1:NDp,1:LNDp,k-1) * F_GE(1:LNDp,k-1) ) ;
            % 
        end
        k = TotNumDM ;
        NDp=DegDM(k)+1; LNDp=DegDM(k-1)+1;
        tmp = inv( BlockL(1:NDp,1:NDp,k) ...
            - SubMat(1:NDp,1:LNDp,k-1) * SupMat_GE(1:LNDp,1:NDp,k-1) );
        
        BlockLInv(1:NDp,1:NDp,k) = tmp;
        
%        F_GE(1:NDp,k) = tmp * ( F(1:NDp,k) ...
%            - SubMat(1:NDp,1:LNDp,k-1) * F_GE(1:LNDp,k-1) ) ;
        
        % backword substitution
%        k = TotNumDM;   ND=DegDM(k);    NDp = ND+1;
%        v(1:NDp,k) = F_GE(1:NDp,k);
%        for k = TotNumDM-1 : -1 : 1
%            NDp=DegDM(k)+1; RNDp=DegDM(k+1)+1;
%            v(1:NDp,k) = F_GE(1:NDp,k) - SupMat_GE(1:NDp,1:RNDp,k) * v(1:RNDp,k+1);
%        end
        
end  % switch

end