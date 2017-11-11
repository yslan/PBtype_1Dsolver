%  Use geometric multigrid to solve 1D (or 2D) Poisson equation


z = 0:N; z=z'/N; % Points on [0,1]


dx = z(2)-z(1);

n  = N-1;
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);

a = 0.5;
P = 0*A; i=0;
for j=1:(n/2);
    i=i+1; P(i  ,j) = 0.5;
    i=i+1; P(i  ,j) = 1.0;
           P(i+1,j) = 0.5;
end;
i=i+1; P(i,j) = 0.5;
P=P(1:i,1:j);
P

