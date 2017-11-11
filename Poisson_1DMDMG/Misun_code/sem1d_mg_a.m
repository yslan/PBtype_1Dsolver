%  Geometric multigrid to solve 1D Poisson equation with SEM
format compact; format shorte


N = 12;  % Polynomial order
E = 30;  % Number of elements in x

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

nbar=size(Abar,1);               %% Dirichlet conditions at each end
R=speye(nbar); R=R(2:(end-1),:); %% R is the BC restriction matrix - Dirichlet
R=speye(nbar); R=R(2:end,:);     %% R is the BC restriction matrix - Neumann
n=size(R,1);

A = R*Abar*R'; B = R*Bbar*R'; %% Dirichlet operators

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
%  Here we set up a multigrid iteration.
%
%  First  pass: Jacobi smoothing
%  Second pass: Overlapping-Schwarz smoothing
%

close

b=r; x0=5.e-1*rand(n,1); 

for kk=1:40;
   r = R*Q'*BL*fL;  %% RHS

   M=diag(sparse(A));  M=inv(diag(M));  % Jacobi smoothing
   sigma = 2/3;
   sigma = input('Input sigma:');

   x=x0; r=b-A*x;

   m_smooth=10;
   for k=1:m_smooth;

       x = x + sigma*M*r;
       r = b-A*x;

   end;

   %%% check error behavior for model problem
       uL = Q*(R'*x);
       eL = sin(pi*xL/2)-uL;
       plot(xL,eL,'r.-')
   %%% 
end;



