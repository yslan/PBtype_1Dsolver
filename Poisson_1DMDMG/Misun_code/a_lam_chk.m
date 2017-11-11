%  Geometric multigrid to solve 1D Poisson equation


%  Usage:   N=40; m=1; mgrid
%
%  N = number of points-1, including domain endpoints.
%      Number of degrees of freedom is N-1.
%
%  k = number of fine-grid smoothings per V-cycle
%

hold off; format compact;  close all


%  Illustration of Damped-Jacobi smoothing effect on error

z = 0:N; z=z'/N; h = z(2)-z(1); % Points and gridspacing on [0,1]
n = N-1; e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n) / (h*h);

[V,L]=eig(full(A)); L=diag(L);
lmin = min(L);
emin = pi*pi-lmin;
e2=12*emin/( (h*h)*(pi^4) );

lmax = max(L);
amax = 4/(h*h);
amax-lmax
emax = ( (4/(h*h))-lmax ) / (pi*pi);

%emax = .5*(h*h) * ( (4/(h*h))-lmax ) ;
%em   = 2*emax/( (h*h)*(pi^2) )







