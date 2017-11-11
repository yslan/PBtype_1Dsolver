function y = MLx(x,DegDM,BlockL,SupMat,SubMat)
%global DegDM
%global BlockL
%global SupMat
%global SubMat
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

TotNumDM = length(DegDM);
ymdarray = zeros(size(BlockL,1),TotNumDM);
ylength = sum(DegDM)+TotNumDM;
y=zeros(ylength,1);
xmdarray = ymdarray;
xmdarray = vec2mdarray(DegDM,x,ymdarray);

switch TotNumDM
    case 1
       k=1; ND=DegDM(k); NDp=ND+1; 
       ymdarray(1:NDp,k) = BlockL(1:NDp,1: NDp) * xmdarray(1:NDp,k);
       y(1:NDp) = ymdarray(1:NDp,k);
 
    otherwise
       k=1; ND=DegDM(k); NDp=ND+1; RND=DegDM(k+1); RNDp = RND+1;
       ymdarray(1:NDp,k) = BlockL(1:NDp,1: NDp) * xmdarray(1:NDp,k)...
                         + SupMat(1:NDp,1:RNDp) * xmdarray(1:RNDp,k+1);
                    
       for k=2:TotNumDM-1
           NDp=DegDM(k)+1; RNDp=DegDM(k+1)+1; LNDp=DegDM(k-1)+1;
           ymdarray(1:NDp,k) = SubMat(1:NDp,1:LNDp,k-1) * xmdarray(1:LNDp,k-1)...
                             + BlockL(1:NDp,1: NDp,k  ) * xmdarray(1: NDp,k  )...
                             + SupMat(1:NDp,1:RNDp,k  ) * xmdarray(1:RNDp,k+1);
       end
       
       k = TotNumDM; NDp = DegDM(k)+1; LNDp = DegDM(k-1)+1; 
       ymdarray(1:NDp,k) = SubMat(1:NDp,1:LNDp,k-1) * xmdarray(1:LNDp,k-1)...
                         + BlockL(1:NDp,1: NDp,k  ) * xmdarray(1: NDp,k  );
       
       y = mdarray2vec(DegDM,ymdarray);
end


end

