
tic
clear all
for N=5:100;

[Ah,Bh,Ch,Dh,z,w] =  semhat(N);

x=z; y=z; [X,Y]=ndgrid(x,y);

n = N-1;  N1=N+1;
P = zeros(N1,n); P(2:(N1-1),:)=speye(n); R=P';

A=R*Ah*R'; B=R*Bh*R';

[S,L] = eig(A,B);

scale = diag(S'*B*S); scale=diag(1./sqrt(scale));
S = S*scale;

sx = sin(pi*X); sy=sin(pi*Y);
cx = cos(pi*X); cy=cos(pi*Y);

g = sx.*sy;
eg = exp(g);
gx=pi*cx.*sy; gy=pi*sx.*cy;
gxx=-pi*pi*sx.*sy; gyy=-pi*pi*sx.*sy;

f=-(gxx+gyy+gx.*gx + gy.*gy).*eg;
uex = eg-1.;

r = R*(Bh*f*Bh')*R';

L=sparse(L); n=size(L,1); I=speye(n);

d=kron(I,L)+kron(L,I);

d=diag(d); d=reshape(1./d,n,n);

u =           S*(d.*(S'*r*S))*S';
u =  R'*u*R;

mesh(X,Y,u)

err = u-uex;
%mesh(X,Y,1.e5*err)

eN(N)=max(max(err));
NN(N)=N;

end;

toc
semilogy(NN,eN,'ro-')


