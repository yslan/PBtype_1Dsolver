%  Geometric multigrid to solve 1D Poisson equation with SEM
format compact; format shorte


N = 4;  % Polynomial order
E = 2;  % Number of elements in x

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

nbar=size(Abar,1); n=nbar-2;     %% Dirichlet conditions at each end
R=speye(nbar); R=R(2:(end-1),:); %% R is the BC restriction matrix - Dirichlet
R=speye(nbar); R=R(2:end,:);     %% R is the BC restriction matrix - Neumann

A = R*Abar*R'; B = R*Bbar*R'; %% Dirichlet operators

%
%  Test 1D code here:  u=sin(pi*x) --> f=pi*pi*sin(pi*x);
%
   fL=pi*pi*sin(pi*xL);             %% Dirichlet
   fL=pi*pi*sin(pi*xL/2)/4;         %% Neumann
   u=A\(R*(Q'*(BL*fL)));
   uL=Q*(R'*u);
   plot(xL,uL,'r.-');
   error=max(abs(uL-sin(pi*xL)));   %% Dirichlet
   error=max(abs(uL-sin(pi*xL/2))); %% Neumann
   [N E error]


%
%
%






