function y = bbC(x,DegDM,Mass)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

TotNumDM = length(DegDM);
ymdarray = zeros(size(Mass,1),TotNumDM);
ylength = sum(DegDM)+TotNumDM;
y=zeros(ylength,1);
xmdarray = ymdarray;
xmdarray = vec2mdarray(DegDM,x,ymdarray);


for k=1:TotNumDM
    NDp=DegDM(k)+1;
    ymdarray(1:NDp,k) = 1./ Mass(1:NDp,k) .* xmdarray(1:NDp,k);
end

y = mdarray2vec(DegDM,ymdarray);

end