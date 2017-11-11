function B_inv = inv_block_mat(A_inv,w,z,loc)
%  Solving the inverse of the matrix 
%   case:loc = 1
%           A   |w
%  B =     _____|_
%           w^T |z
%  B: (n+1)*(n+1)
%  A: n*n, 
%  w: n*1, 
%  z: 1*1
%
%   case: loc = 2
%           z |   w^T
%  B =     ---|----------
%           w |     A 
%             |    
%  
%  B: (n+1)*(n+1)
%  A: n*n, 
%  w: n*1, 
%  z: 1*1 

[N,~] = size(A_inv);

B_inv = zeros(N+1,N+1);
k = z - w'*A_inv*w;
k_inv = 1/k;
switch loc
    case 1
%         B(1:N,1:N) = A_inv^(-1);
%         B(N+1,1:N) = w';
%         B(1:N,N+1) = w;
%         B(N+1,N+1) = z;
        
        B_inv(1:N,1:N) =  A_inv + k_inv *(A_inv*w) * (w'*A_inv);
        B_inv(N+1,1:N)   = -k_inv * w'*A_inv;
        B_inv(1:N,N+1)   = -k_inv * A_inv*w;
        B_inv(N+1,N+1) =  k_inv;
    case 2
%         B(1,1) = z;
%         B(1,2:N+1) = w';
%         B(2:N+1,1) = w;
%         B(2:N+1,2:N+1) =  A_inv^(-1);
        
        B_inv(2:N+1,2:N+1) =  A_inv + k_inv *(A_inv*w) * (w'*A_inv);
        B_inv(1,2:N+1)     = -k_inv * w'*A_inv;
        B_inv(2:N+1,1)     = -k_inv * A_inv*w;
        B_inv(1,1)         =  k_inv;
end











