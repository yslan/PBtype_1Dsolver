function mdarray = vec2mdarray(DegDM,vec,ref_mdarray)

TotNumDM = length(DegDM);
NDmax = max(DegDM); NDmaxp = NDmax+1;

mdarray = zeros(size(ref_mdarray));

k = 1; ND = DegDM(k); NDp = ND+1;
mdarray(1:NDp,1) = vec(1:NDp);

for k=2:TotNumDM
    
    ND = DegDM(k); NDp = ND+1; 
    StartIndex = sum(DegDM(1:k-1))+(k-1)+1;
    EndIndex = sum(DegDM(1:k)) + k;

    mdarray(1:NDp,k) = vec(StartIndex:EndIndex);

end