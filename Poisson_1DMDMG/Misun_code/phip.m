      function[p] =  phip(x,N) % Compute bubble functions up to degree N

      [p]=legendre(x,N);

      for k=N:-1:3; i=k+1;
          p(:,i) = p(:,i)-p(:,i-2);
      end;
