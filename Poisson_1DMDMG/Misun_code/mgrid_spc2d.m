
n=50; r=1;

dx = 1/(n+1);

k=1:n; k=k'; l=k;

[K,L]=ndgrid(k,l);

lkl = K;
for i=1:n; for j=1:n;
 lkl(i,j)=(2/(dx*dx))*( (1-cos(i*pi*dx)) + (1-cos(j*pi*dx))/(r*r));
end; end;

mesh(K,L,lkl)





