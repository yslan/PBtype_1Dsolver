      function[Ah,Bh,Ch,Dh,w] =  femhatg(z)
%
%     1D fem Stiffness, Mass, and Convection matrix analogies to SEMhat
%

      N  = size(z,1)-1;
      h  = diff(z);

      n  = N+1;

      Ep = diag(ones(n-1,1),1);  Ep = sparse(Ep); 
      Em = Ep';
      D  = speye(n,n);
      Ah = 0*(Ep+Em+D);   % Tridiagonal Ah placeholder
      Bh = Ah;        % Tridiagonal Bh placeholder
      Dh = Ah;        % Tridiagonal Bh placeholder
      Ch = Ah;        % Tridiagonal Bh placeholder

      for e=1:N; 
        Ah(e+0,e+0) = Ah(e+0,e+0) + 1./h(e);
        Ah(e+0,e+1) = Ah(e+0,e+1) - 1./h(e);
        Ah(e+1,e+0) = Ah(e+1,e+0) - 1./h(e);
        Ah(e+1,e+1) = Ah(e+1,e+1) + 1./h(e);

        Bh(e+0,e+0) = Bh(e+0,e+0) + h(e)/6;
        Bh(e+0,e+1) = Bh(e+0,e+1) + h(e)/6;
        Bh(e+1,e+0) = Bh(e+1,e+0) + h(e)/6;
        Bh(e+1,e+1) = Bh(e+1,e+1) + h(e)/6;

      end;

      w     = Bh*ones(n,1);

