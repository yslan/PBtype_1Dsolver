%  Use geometric multigrid to solve 1D (or 2D) Poisson equation


z = 0:N; z=z'/N; % Points on [0,1]


dx = z(2)-z(1);

n  = N-1;
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);




ub = 0*z;  ue=[0; ue; 0];  % Extended solution

%%% Smoother
S = diag(A); S=(1.+0*S)./S; S=diag(S); sigma = 0.6666;

ue = rand(n,1);  % Exact solution (Endpoints assumed to be 0)
b = A*ue;        % RHS
x = 0*b;

for k=1:40;

   r = b-A*x;
   for j=1:3

%      Diagnositics
       xe = [0; x; 0];
       ee = ue - xe;
       plot(z,xe,'b-',z,ee,'r-',z,ue,'k--')
       pause(0.1)

       dx = sigma*(S*r);
       x  = x + dx;
       r  = r - A*dx;
   end;

   pause






end;





























uuuu
%  Use geometric multigrid to solve 1D (or 2D) Poisson equation


z = 0:N; z=z'/N; % Points on [0,1]


dx = z(2)-z(1);

n  = N-1;
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);

P = 0*A; i=0;
for j=1:(n/2);
    i=i+1; P(i  ,j) = 0.5; i=i+1; P(i  ,j) = 1.0; P(i+1,j) = 0.5;
end;
i=i+1; P(i,j) = 0.5; P=P(1:i,1:j);

