
clear all; close all

tic
for N=4:280;

[Ah,Bh,Ch,Dh,z,w] =  semhat(N);

x=z; y=z; [X,Y]=ndgrid(x,y);

n = N-1;  N1=N+1;
P = zeros(N1,n); P(2:(N1-1),:)=speye(n); R=P';

A=R*Ah*R'; B=R*Bh*R';
[S,L] = eig(A,B); n=size(S,1); 

scale = diag(S'*B*S); scale=diag(1./sqrt(scale)); S = S*scale;

I=speye(n); L=sparse(L);

d = kron(I,L)+kron(L,I); 
d=diag(d); d=1./d;
d = reshape(d,n,n);

sx=sin(pi*X); sy=sin(pi*Y); cx=cos(pi*X); cy=cos(pi*Y);
g=sx.*sy; 
gxx = -pi*pi*g; gyy = -pi*pi*g; 
gx=pi*cx.*sy; gy=pi*sx.*cy;

eg  = exp(g);
uex = eg-1.0;
f=-(gxx+gyy+gx.*gx+gy.*gy).*eg;

r = (R*Bh)*f*(R*Bh)';

u = S*(d.*(S'*r*S))*S';
u = R'*u*R;

%  mesh(X,Y,uex); pause
%  mesh(X,Y,u); pause

err = abs(u-uex);
% mesh(X,Y,err); drawnow; pause(1)
emx = max(max(err));
NN(N)=N;
eN(N)=emx;
end;

semilogy(NN,eN,'ro-')

toc
