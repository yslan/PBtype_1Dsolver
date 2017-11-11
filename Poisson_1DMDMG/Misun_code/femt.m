
[Ah,Bh,Ch,Dh,z,w]=semhat(N);
[Af,Bf,Cf,Df,w]=femhatg(z); 

n=size(Af,1);
P=eye(n-2); P=[ 0*P(1,:) ; P ; 0*P(1,:) ]; R=P';

A=R*Af*R'; B=R*Bf*R';


%[V,L]=eig(full(A),full(B));
[V,L]=eig(full(A));

L=diag(L); 
min(L)/(pi*pi/4)

kappa1=kappa;
kappa=max(L)/min(L)

ratio = kappa1/kappa

Ah=R*Ah*R'; Bh=R*Bh*R';
[V,L]=eig(full(Ah),full(A)); L=diag(L);
kapp2=(max(L)/min(L))/(pi*pi/4)


