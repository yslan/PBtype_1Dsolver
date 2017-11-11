%  Geometric multigrid to solve 1D Poisson equation with SEM
format compact; format shorte


N = 128;  % Polynomial order
E = 1;  % Number of elements in x

Lx = 1.0; % Domain size

Le=ones(E,1)*(Lx/E);  % Uniform element size

[Ah,Bh,Ch,Dh,z,w]=semhat(N);
[Q]=semq(E,N,0);

% CONSTUCT GLOBAL OPERATORS

Ie=speye(E);

AL=diag(2./Le);AL=kron(AL,Ah); Abar=Q'*AL*Q;   % Abar = Neumann operator
BL=diag(Le/2); BL=kron(BL,Bh); Bbar=Q'*BL*Q; 

xL = zeros(N+1,E); x0=0;
for e=1:E; xL(:,e)=x0 + Le(e)*(z+1)/2; x0=xL((N+1),e); end; 
xL=reshape(xL,(N+1)*E,1);


% Coarse grid setting
Nc=floor(N/2); %% Piecewise linear coarse-grid space
Nc=floor(2)  ; %% Piecewise linear coarse-grid space
[Ac,Bc,Cc,Dc,zc,wc]=semhat(Nc);
Jh = interp_mat(z,zc);
[Qc]=semq(E,Nc,0);
JL = kron(Ie,Jh);
mult = 1+0*xL; mult=Q*Q'*mult; multL=1./mult;
J=zeros(N*E+1,Nc*E+1); J=sparse(J);
i0=1;j0=1;
for e=1:E;
  i1=i0+N; j1=j0+Nc;
  J(i0:i1,j0:j1)=Jh;
  i0=i1; j0=j1;
end;
Acbar=J'*Abar*J;
Bcbar=J'*Bbar*J;


nbar=size(Abar,1);
R=speye(nbar); R=R(2:(end-1),:); %% R is the BC restriction matrix - Dirichlet
R=speye(nbar); R=R(2:end,:);     %% R is the BC restriction matrix - Neumann
n=size(R,1);

ncbar=size(Acbar,1);
Rc=speye(ncbar); Rc=Rc(2:(end-1),:); %% Rc is the BC restriction matrix - Dirichlet
Rc=speye(ncbar); Rc=Rc(2:end,:);     %% Rc is the BC restriction matrix - Neumann
nc=size(Rc,1);

A  = R*Abar*R'; B = R*Bbar*R';        %% Dirichlet operators
Ac = Rc*Acbar*Rc'; Bc = Rc*Bcbar*Rc'; %% Dirichlet operators
J  = R*J*Rc';
%
%  Test 1D code here:  u=sin(pi*x) --> f=pi*pi*sin(pi*x);
%
   fL=pi*pi*sin(pi*xL);             %% Dirichlet
   fL=pi*pi*sin(pi*xL/2)/4;         %% Neumann
   r =R*Q'*BL*fL;                   %% RHS
   u=A\r;
   uL=Q*(R'*u);
   plot(xL,uL,'r.-');
   error=max(abs(uL-sin(pi*xL)));   %% Dirichlet
   error=max(abs(uL-sin(pi*xL/2))); %% Neumann
   [N E error]

%
%  Here, set up the restriction matrices for overlapping Schwarz
%
i0=1; No=N+1;
for e=1:E;
%   i1=e*(N+2); i0=i1-(N+3); 
    i1=e*(N+3); i0=i1-(N+4); 
    ischwarz(1,e)=max(1,i0);
    ischwarz(2,e)=min(n,i1);
end;

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
%  sigma = input('Input sigma:');

   x=x0; r=b-A*x;

   m_smooth=2;
   for k=1:m_smooth;

%      x = x + sigma*M*r;

       z = 0*x;
       for e=1:E;
         i0=ischwarz(1,e); i1=ischwarz(2,e);
         z(i0:i1)=z(i0:i1)+A(i0:i1,i0:i1)\r(i0:i1);
       end;

       x = x + sigma*z;
       r = b-A*x;
   end;

   %%% check error behavior for model problem
       uL = Q*(R'*x);
       eL = sin(pi*xL/2)-uL;
       plot(xL,eL,'r.-')
   %%% 

%% Coarse-grid correction  x <--- x + e, where e is approximated on coarse grid.

   rc = J'*r;    % Restrict residual to coarse grid
   ec = Ac\rc;   % Solve coarse-grid problem:  Ac*ec = rc;
   ef = J*ec;    % Interpolate coarse-grid error to fine grid

   x = x+ef;

   %%% check error behavior for model problem
       uL = Q*(R'*x);
       eL = sin(pi*xL/2)-uL;
       hold on; plot(xL,eL,'b-'); hold off
       pause
   %%% 
   x0=x;

end;

