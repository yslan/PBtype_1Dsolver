function  v = MDGaussSol2(DegDM,F,BlockLInv,SupMat_GE,SubMat)
%Gauss Elimination 
%Forward sweeping

TotNumDM = length(DegDM); 
switch TotNumDM
    case 1
        k = 1; ND = DegDM(k); NDp = ND+1;
        v(1:NDp,k) = BlockLInv(1:NDp,1:NDp,k) * F(1:NDp,k);
        
    otherwise
        % forward sweeping
        k = 1; NDp=DegDM(k)+1; RNDp=DegDM(k+1)+1;
        F_GE(1:NDp,k) = BlockLInv(1:NDp,1:NDp,k) * F(1:NDp,k) ; 
        
        for k = 2 : TotNumDM-1
            LNDp = DegDM(k-1)+1; NDp = DegDM( k )+1; RNDp = DegDM(k+1)+1;
            
            F_GE(1:NDp,k) = BlockLInv(1:NDp,1:NDp,k) * ( F(1:NDp,k) ...
                - SubMat(1:NDp,1:LNDp,k-1) * F_GE(1:LNDp,k-1) ) ;
            
        end
        k = TotNumDM ;
        NDp=DegDM(k)+1; LNDp=DegDM(k-1)+1;
        
        F_GE(1:NDp,k) = BlockLInv(1:NDp,1:NDp,k) * ( F(1:NDp,k) ...
            - SubMat(1:NDp,1:LNDp,k-1) * F_GE(1:LNDp,k-1) ) ;
        
        % backward sweeping
        k = TotNumDM;   ND=DegDM(k);    NDp = ND+1;
        v(1:NDp,k) = F_GE(1:NDp,k);
        for k = TotNumDM-1 : -1 : 1
            NDp=DegDM(k)+1; RNDp=DegDM(k+1)+1;
            v(1:NDp,k) = F_GE(1:NDp,k) - SupMat_GE(1:NDp,1:RNDp,k) * v(1:RNDp,k+1);
        end
        
end  % switch

end