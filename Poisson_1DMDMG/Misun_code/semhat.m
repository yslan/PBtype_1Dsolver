      function[Ah,Bh,Ch,Dh,z,w] =  semhat(N)
%
%                                                 ^
%     Compute the single element 1D SEM Stiffness Mass, and Convection
%     matrices, as well as the points and weights for a polynomial
%     of degree N
%

      [z,w] = zwgll(N);

      Bh    = diag(w);
      Dh    = dhat(z);

      Ah    = Dh'*Bh*Dh;
      Ch    = Bh*Dh;

