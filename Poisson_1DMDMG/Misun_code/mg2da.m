%  Geometric multigrid to solve 2D Poisson equation

%  Usage:   N=40; k=1; mg2d
%
%  N = number of points-1, including domain endpoints.
%      Number of degrees of freedom is N-1.
%
%  k = number of fine-grid smoothings per V-cycle
%

msmooth = k;

z  = 0:N; z=z'/N;  % Points on [0,1]

lenx=1; leny=5.0;
x = z*lenx; dx = x(2)-x(1); 
y = z*leny; dy = y(2)-y(1); 

[X,Y]=ndgrid(x,y);
[nx,ny]=size(X); nx=nx-2; ny=ny-2; % ASSUMES DIRICHLET CONDITIONS in x and y!

n  = N-1;
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);



%  Illustration of Damped-Jacobi smoothing effect on error

z  = 0:N; z=z'/N; h = z(2)-z(1); % Points and gridspacing on [0,1]
n  = N-1; e = ones(n,1);
Ax = spdiags([-e 2*e -e], -1:1, n, n) / (dx*dx); Bx = speye(n); Dx=diag(diag(Ax));
Ay = spdiags([-e 2*e -e], -1:1, n, n) / (dy*dy); By = speye(n); Dy=diag(diag(Ay));
A = kron(By,Ax)+kron(Ay,Bx);
D = reshape(diag(A),nx,ny); sigma=0.66666;    % Smoother

[Vx,Lx]=eig(full(Dx\Ax)); Lx=eig(Lx);  %% Note: eigenfunctions indeterminate up to a sign
Sx=Vx'*Vx; Sx=diag(Sx); Sx=1./sqrt(Sx); Sx=diag(Sx); Vx=Vx*Sx; 
for k=1:n; if Vx(1,k) < 0; Vx(:,k)=-Vx(:,k); end; end;

[Vy,Ly]=eig(full(Dy\Ay)); Ly=eig(Ly);  %% Note: eigenfunctions indeterminate up to a sign
Sy=Vy'*Vy; Sy=diag(Sy); Sy=1./sqrt(Sy); Sy=diag(Sy); Vy=Vy*Sy; 
for k=1:n; if Vy(1,k) < 0; Vy(:,k)=-Vy(:,k); end; end;


ex= rand(nx,ny);
ex= 1+0*ex; ex=Vx*ex*Vy';
f = Ax*ex + ex*Ay';            % f = A*exact
u = 0*f; r = f-(Ax*u+u*Ay');   % Initial guess and residual
mesh(ex); pause

v  = 0; ispec=1; % ispec=1 --> plot spectra
ip = plot_err2d(u,ex,X,Y,Vx,Vy,Lx,Ly,'g-',v,0,ispec); pause

m=2;
for v=1:10;                    % 10 V-cycles
  for j=1:msmooth;             % Apply smoother m times
    u=u+sigma*(D.\r); r=f-(Ax*u+u*Ay');
    hold off
    ip = plot_err2d(u,ex,X,Y,Vx,Vy,Lx,Ly,'b-',j,0,ispec); 
    pause
  end;
  ip = plot_err2d(u,ex,X,Y,Vx,Vy,Lx,Ly,'b-',v,0,ispec); 

  ef = R'*(Ac \ (R*r) );       % Coarse-grid correction
  u  = u + ef;
  ip = plot_err2d(u,ex,X,Y,Vx,Vy,Lx,Ly,'r-',v,0,ispec);  pause

end;

