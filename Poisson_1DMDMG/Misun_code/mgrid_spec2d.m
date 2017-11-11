
n=50; r=1;

dx = 1/(n+1);

k=1:n; k=k'; l=k;

[K,L]=ndgrid(k,l);

lkl = K;
for i=1:n; for j=1:n;
 lkl(i,j)=(2/(dx*dx))*( (1-cos(i*pi*dx)) + (1-cos(j*pi*dx))/(r*r));
end; end;

mesh(K,L,lkl)


%  Geometric multigrid to solve 1D Poisson equation


%  Usage:   N=40; m=1; mgrid
%
%  N = number of points-1, including domain endpoints.
%      Number of degrees of freedom is N-1.
%
%  k = number of fine-grid smoothings per V-cycle
%

hold off; format compact;  close all

z  = 0:N; z=z'/N;  % Points on [0,1]


dx = z(2)-z(1);

n  = N-1;
e = ones(n,1);
Ax = spdiags([-e 2*e -e], -1:1, n, n);
Ix = speye(n);

r = 10;
Ay= Ax/(r*r); Iy=Ix;

A=kron(Iy,Ax)+kron(Ay,Ix);
D=diag(diag(A));

[V,L]=eig(full(D\A)); L=eig(L);  %% Note: eigenfunctions indeterminate up to a sign
S=V'*V; S=diag(S); S=1./sqrt(S); S=diag(S); V=V*S; 




