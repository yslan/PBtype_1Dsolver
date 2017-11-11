%  Geometric multigrid to solve 1D Poisson equation

%  Usage:   N=40; k=1; mgrid
%
%  N = number of points-1, including domain endpoints.
%      Number of degrees of freedom is N-1.
%
%  k = number of fine-grid smoothings per V-cycle
%

z  = 0:N; z=z'/N;  % Points on [0,1]


dx = z(2)-z(1);

n  = N-1;
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);


%  Illustration of Damped-Jacobi smoothing effect on error

z = 0:N; z=z'/N; h = z(2)-z(1); % Points and gridspacing on [0,1]
n = N-1; e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n) / (h*h);
D = diag(diag(A)); sigma=0.66666;    % Smoother


Nc = N/2; zc = 0:Nc; zc=zc'/Nc;     %% Build Level 2 Operators:
P = 0*A; i=0;                       %% 
for j=1:(n/2); i=i+1; P(i,j)=.5; i=i+1; P(i,j)=1; P(i+1,j)=.5; end;
               i=i+1; P(i,j)=.5;
P=P(1:i,1:j);  %% PROLONGATION OPERATOR
R=P';          %% RESTRICTION OPERATOR
Ac = R*A*R';   %% COARSE-GRID SYSTEM

ex= rand(n,1); f = A*ex;   % f = A*exact
ex= [0; ex; 0];  % Extend exact solution to boundaries (for plotting)

u = 0*f; r = f-A*u;            % Initial guess and residual

for V=1:10;                    % 10 V-cycles
  for j=1:k                    % Apply smoother k times
    u=u+sigma*(D\r); r=f-A*u;
  end;
  ef = R'*(Ac \ (R*r) );       % Coarse-grid correction
  u  = u + ef;

  %DIAGNOSITICS: Use known solution to show error behavior.
  ue=[0;u;0];ee=ex-ue;plot(z,ue,'bo-',z,ee,'r-',z,ex,'k--',z,0*z,'k-'); 
  strg=sprintf('%s','Error after ',int2str(V),' V-cycles');
  xlabel('x'); ylabel('Solution / Error'); title(strg); axis square
  legend('numeric','error','exact'); set (gcf,'color',[1 1 1]);
  norm(ee), pause
end;
