
clear all; close all; format compact

for k=1:22; scale = k*(k-1);

 for N=3:40;
   [Ah,Bh,Ch,Dh,z,w] = semhat(N);
   n=size(Ah,1); R=eye(n); R=R(2:n-1,:); %Restriction matrix for Dirichlet BCs

   A=R*Ah*R';  Bb=R*Bh;

   if k<3;              %% An analytic function
      ue=cos(.5*pi*z);
      f =0.25*pi*pi*cos(0.5*pi*z);
   else;
      ue=1-abs(z).^k;   %% A function with limited regularity.
      f =scale*abs(z).^(k-2);
   end;

   u =R'*(A\(Bb*f));
   err=u-ue;
   e = max(abs(err));

   en(N)=e;
   en(N)=sqrt(err'*Bh*err)/2;
   nn(N)=N;

 %  [N e]
 %  plot(z,ue,'kx-',z,u,'ro-'); pause(.1)
 end;

 k
 if k < 3;
   if mod(k,2)==0; figure(1); loglog(nn,en,'b.-'); hold on; pause(.1); else;
                       figure(2); loglog(nn,en,'b.-'); hold on; pause(.1); end;
 else;
   if mod(k,2)==0; figure(1); loglog(nn,en,'r.-'); hold on; pause(.1); else;
                       figure(2); loglog(nn,en,'k.-'); hold on; pause(.1); end;
 end;

 if k==5; pause; end;

end;

