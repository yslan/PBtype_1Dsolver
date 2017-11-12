function [f1,f2] = userf(C1,C2)
steg = [1,1.2;1.2,1];

f1 = steg(1,1)*C1 + steg(1,2)*C2;
f2 = steg(2,1)*C1 + steg(2,2)*C2;

% f1 = zeros(size(c1));
% f2 = zeros(size(c2));

end
