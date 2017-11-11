     function [A,B,C,X] = sem1c(E,N,Lx,bc);  
%
%    Boundary conditions - 0 = Dirichlet
%                          2 = Periodic
%
%    Build SEM A & B matrix for the 1D Laplacian with
%    Dirichlet or boundary conditions on the interval [0:Lx].
%
%    E = number of spectral elements
%    N = polynomial order of spectral elements
%
%    A = disrete Laplacian  (SPD, null space = constant mode)
%    C = disrete convection operator
%    B = mass matrix
%    X = nodal spacing
%
%    Ah = disrete Laplacian for single element
%    Bh = mass matrix for single element
%    Ch = convection operator for single element
%    Dh = derivative matrix for single element
%
%

 [Ah,Bh,Ch,Dh,z,w] = semhat(N);

%
%
%   Assemble local operators
%
%
 Ah = sparse(Ah);
 Bh = sparse(Bh);
%
 Ah = sparse(Ah);
 Bh = sparse(Bh);
 Ch = sparse(Bh*Dh);

 Lxe   = Lx/E;

 if bc==2;
    Ie = sparse(eye(E));
    Ue = sparse(diag(ones(E-1,1),1));
    Ue(E,1) = 1;
    Le = Ue';

    n1      = N+1;

    Ad      = Ah(1:N,1:N); Ad(1,1) = Ad(1,1) + Ah(n1,n1);
    Au      = zeros(N,N);  Au(:,1) = Ah(1:N,n1);
    Al      = zeros(N,N);  Al(1,:) = Ah(n1,1:N);
    A       = (2./Lxe)*(kron(Ie,Ad) + kron(Ue,Au) + kron(Le,Al));

    Cd      = Ch(1:N,1:N); Cd(1,1) = Cd(1,1) + Ch(n1,n1);
    Cu      = zeros(N,N);  Cu(:,1) = Ch(1:N,n1);
    Cl      = zeros(N,N);  Cl(1,:) = Ch(n1,1:N);
    C       = (kron(Ie,Cd) + kron(Ue,Cu) + kron(Le,Cl));

    Bd      = Bh(1:N,1:N); Bd(1,1) = Bd(1,1) + Bh(n1,n1);
    B       = (Lxe/2.)*kron(Ie,Bd);
    n       = size(B,1);
    X       = zeros(n,1);
    for e=1:E;
        ii = 1 + N*(e-1);
        X(ii:ii+N)         = Lxe*(e-1) + (Lxe/2)*(z+1);
    end;
 end;

 if bc < 2;
    Ah = (2/Lxe)*Ah;
    Bh = (Lxe/2)*Bh;
    Ch = Bh*Dh;
    n  = E*N+1;

    A  = 0.*speye(n); B=A; C=A;
    X  = zeros(n,1);
    for e=1:E;
        ii = 1 + N*(e-1);
        A(ii:ii+N,ii:ii+N) = A(ii:ii+N,ii:ii+N) + Ah;
        B(ii:ii+N,ii:ii+N) = B(ii:ii+N,ii:ii+N) + Bh;
        C(ii:ii+N,ii:ii+N) = C(ii:ii+N,ii:ii+N) + Ch;
        X(ii:ii+N)         = Lxe*(e-1) + (Lxe/2)*(z+1);
    end;
 %  spy(A); pause
    if bc == 0;
       A = A(2:n-1,2:n-1);
       B = B(2:n-1,2:n-1);
       C = C(2:n-1,2:n-1);
       X = X(2:n-1);
    end;
 end;
%spy(A); pause
