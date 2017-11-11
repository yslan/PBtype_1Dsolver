function  [v]=GaussEliminationMD(TotNumDM,DegDM,F,BlockL,SubMat,SupMat)
%Gauss Elimination 
%Forward sweeping

switch TotNumDM
    case 1
        k = 1; ND = DegDM(k); NDp = ND+1;
        tmp = inv( BlockL(1:NDp,1:NDp,k)) ;
        v(1:NDp,k) = tmp * F(1:NDp,k) ;
        
    otherwise
        
        k = 1; NDp=DegDM(k)+1; RNDp=DegDM(k+1)+1;
        tmp = inv( BlockL(1:NDp,1:NDp,k)) ;
        SupMat_GE(1:NDp,1:RNDp,k) = tmp * SupMat(1:NDp,1:RNDp,k) ;
        F_GE(1:NDp,k) = tmp * F(1:NDp,k) ;
        
        for k = 2 : TotNumDM-1
            LNDp = DegDM(k-1)+1; NDp = DegDM( k )+1; RNDp = DegDM(k+1)+1;
            
            tmp = inv( BlockL(1:NDp,1:NDp,k) ...
                - SubMat(1:NDp,1:LNDp,k-1) * SupMat_GE(1:LNDp,1:NDp,k-1) );
            
            SupMat_GE(1:NDp,1:RNDp,k) = tmp * SupMat(1:NDp,1:RNDp,k) ;
            
            F_GE(1:NDp,k) = tmp * ( F(1:NDp,k) ...
                - SubMat(1:NDp,1:LNDp,k-1) * F_GE(1:LNDp,k-1) ) ;
            
        end
        k = TotNumDM ;
        NDp=DegDM(k)+1; LNDp=DegDM(k-1)+1;
        tmp = inv( BlockL(1:NDp,1:NDp,k) ...
            - SubMat(1:NDp,1:LNDp,k-1) * SupMat_GE(1:LNDp,1:NDp,k-1) );
        
        F_GE(1:NDp,k) = tmp * ( F(1:NDp,k) ...
            - SubMat(1:NDp,1:LNDp,k-1) * F_GE(1:LNDp,k-1) ) ;
        
        % backword substitution
        k = TotNumDM;   ND=DegDM(k);    NDp = ND+1;
        v(1:NDp,k) = F_GE(1:NDp,k);
        for k = TotNumDM-1 : -1 : 1
            NDp=DegDM(k)+1; RNDp=DegDM(k+1)+1;
            v(1:NDp,k) = F_GE(1:NDp,k) - SupMat_GE(1:NDp,1:RNDp,k) * v(1:RNDp,k+1);
        end
        
end  % switch

end