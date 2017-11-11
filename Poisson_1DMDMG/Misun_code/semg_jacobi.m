%  Geometric multigrid to solve 1D Poisson equation with SEM
clear all
format compact; format shorte


N = 14;  % Polynomial order
E = 13;  % Number of elements in x
N =  7;  % Polynomial order
E =  2;  % Number of elements in x

Lx= 1.0; % Domain size

Le= ones(E,1)*(Lx/E);  % Uniform element size

% Ah: symmetric positive definite, Ch: skew symmetric, Dh: skew-centrosymmetric; 
[Ah,Bh,Ch,Dh,z,w]=semhat(N); 

[Q]=semq(E,N,0);  % check structure of Q with full(Q)

% CONSTUCT GLOBAL OPERATORS

Ie=speye(E);

% Jacobian J= Le/2;
AL=diag(2./Le);AL=kron(AL,Ah); Abar=Q'*AL*Q;   
BL=diag(Le/2); BL=kron(BL,Bh); Bbar=Q'*BL*Q; 

xL = zeros(N+1,E); x0=0;
for e=1:E; xL(:,e)=x0 + Le(e)*(z+1)/2; x0=xL((N+1),e); end; 
xL=reshape(xL,(N+1)*E,1);


Nc=floor(N/2); % Piecewise linear coarse-grid space
[Ac,Bc,Cc,Dc,zc,wc]=semhat(Nc);

%interpolation of fine grid on the coarse grid exapnsion basis on reference domain
Jh = interp_mat(z,zc); 
[Qc]=semq(E,Nc,0);
JL = kron(Ie,Jh); % interpolation on local grids 10 x 6 (N=4, Nc=2)

% multiplicity of grid points
% mult = 1+0*xL; mult=Q*Q'*mult; multL=1./mult;
% interpolation on physical domain using gloabl grids J: 9 x 5 (N=4,Nc=2) 
J=zeros(N*E+1,Nc*E+1); J=sparse(J); 

i0=1;j0=1;
for e=1:E;
  i1=i0+N; j1=j0+Nc;
  J(i0:i1,j0:j1)=Jh;
  i0=i1; j0=j1;
end;
Acbar=J'*Abar*J  ; % Ac on global: (5x9) (9x9) (9x5)  
Bcbar=J'*Bbar*J ;


nbar=size(Abar,1);
R=speye(nbar); R=R(2:(end-1),:); % R is the BC restriction matrix - Dirichlet
R=speye(nbar); R=R(2:end,:);     % R is the BC restriction matrix - Neumann
n=size(R,1);

ncbar=size(Acbar,1);
Rc=speye(ncbar); Rc=Rc(2:(end-1),:); % Rc is the BC restriction matrix - Dirichlet
Rc=speye(ncbar); Rc=Rc(2:end,:);     % Rc is the BC restriction matrix - Neumann
nc=size(Rc,1);

A  = R*Abar*R'; B = R*Bbar*R';        % Dirichlet operators
Ac = Rc*Acbar*Rc'; Bc = Rc*Bcbar*Rc'; % Dirichlet operators
J  = R*J*Rc';

%
%  Test 1D code here:  u=sin(pi*x) --> f=pi*pi*sin(pi*x);
%
   fL=pi*pi*sin(pi*xL);             % Dirichlet
   fL=pi*pi*sin(pi*xL/2)/4;         % Neumann
   r =R*Q'*BL*fL;                   % RHS
   u=A\r;
   uL=Q*(R'*u);
   plot(xL,uL,'r.-');
   error=max(abs(uL-sin(pi*xL)));   % Dirichlet
   error=max(abs(uL-sin(pi*xL/2))); % Neumann
   [N E error]


%
%  Here we set up a multigrid iteration.
%
%  First  pass: Jacobi smoothing
%  Second pass: Overlapping-Schwarz smoothing
%

close

b=r; x0=5.e-1*rand(n,1); 

for vcycle=1:10;

   M=diag(sparse(A));  M=inv(diag(M));  % Jacobi smoothing
   sigma = 2/3;
   sigma =  .3;
%  sigma = input('Input sigma:');

   x=x0; r=b-A*x;

   m_smooth=3;
   for k=1:m_smooth;
       x = x + sigma*M*r;
       r = b-A*x;
   end;

   % check error behavior for model problem
       uL = Q*(R'*x);
       eL = sin(pi*xL/2)-uL;
       plot(xL,eL,'r.-')
   %

% Coarse-grid correction  x <--- x + e, where e is approximated on coarse grid.

   rc = J'*r;    % Restrict residual to coarse grid
   ec = Ac\rc;   % Solve coarse-grid problem:  Ac*ec = rc;
   ef = J*ec;    % Interpolate coarse-grid error to fine grid

   x = x+ef;

   % check error behavior for model problem
       uL = Q*(R'*x);
       eL = sin(pi*xL/2)-uL;
       hold on; plot(xL,eL,'b-'); hold off
       pause
   % 
   x0=x;

end;

