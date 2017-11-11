%  Use geometric multigrid to solve 1D (or 2D) Poisson equation


z = 0:N; z=z'/N; % Points on [0,1]


dx = z(2)-z(1);

n  = N-1;
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);

P = 0*A; i=0;  %% PROLONGATION OPERATOR
for j=1:(n/2);
    i=i+1; P(i  ,j) = 0.5; i=i+1; P(i  ,j) = 1.0; P(i+1,j) = 0.5;
end;
i=i+1; P(i,j) = 0.5;
P=P(1:i,1:j);
R=P';          %% RESTRICTION OPERATOR

Ac = R*A*R';

%%% Smoother
S = diag(A); S=(1.+0*S)./S; S=diag(S); sigma = 0.6666;

ue = 1:n; ue=ue'; ue=(-1).^ue;
ue = 1:n; ue=ue'/(n+1); ue=sin(5*pi*ue);
ue = rand(n,1);  % Exact solution (Endpoints assumed to be 0)
b = A*ue;        % RHS
x = 0*b;

ue = [0;ue; 0];

for k=1:200;

   r = b-A*x;
   for j=1:3

%      Diagnositics
       xe = [0; x; 0]; ee = ue - xe; hold off
       plot(z,xe,'bo-',z,ee,'r-',z,ue,'k--'); 
       legend('numeric','error','exact')
       pause(0.1)

       dx = sigma*(S*r);
       x  = x + dx;
       r  = r - A*dx;   % r = A*error
   end;
   rc = R*r;            % Restrict residual
   ec = Ac\rc;          % Estimate error on coarse grid
   ef = R'*ec;          % Interpolate error estimate to fine grid
   x  = x + ef;         % Add error estimate to current iterate

   efe = [0; ef; 0]; 
   xe  = [0; x; 0]; ee = ue - xe; hold on
   plot(z,xe,'go-',z,efe,'m.-'); 

   pause

end;

uuu
%  Use geometric multigrid to solve 1D (or 2D) Poisson equation


z = 0:N; z=z'/N; % Points on [0,1]


dx = z(2)-z(1);

n  = N-1;
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);

a = 0.5;
