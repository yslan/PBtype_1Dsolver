function [ vec ] = mdarray2vec(DegDM,mdarray)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

TotNumDM = length(DegDM);
veclength = sum(DegDM)+TotNumDM;
vec = zeros(veclength,1);

switch TotNumDM
    case 1
        k=1; ND = DegDM(k); NDp = ND+1;
        vec(1:NDp) = mdarray(1:NDp,k);
        
    otherwise
        
        k = 1; NDp = DegDM(k)+1; 
        vec(1:NDp) = mdarray(1:NDp,1);

        for k=2:TotNumDM
    
            NDp = DegDM(k)+1; 
            StartIndex = sum(DegDM(1:k-1))+(k-1)+1;
            EndIndex = sum(DegDM(1:k)) + k;

            vec(StartIndex:EndIndex) = mdarray(1:NDp,k);
        end
end
        
end

