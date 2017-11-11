function Xlong = StackVec(DegDM,x)

TotNumDM = length(DegDM);
LengthXlong = sum((DegDM+1));
Xlong = zeros(LengthXlong,1);

switch TotNumDM
    case 1
        k=1; ND = DegDM(k); NDp = ND+1;
        Xlong(1:NDp) = x(1:NDp,k);
        
    otherwise
        
        Xlong = [x(1:DegDM(1)+1,1)];
        for k = 2 : TotNumDM 
            Xlong  = [Xlong ; x(1:DegDM(k)+1,k)];
        end
end

end