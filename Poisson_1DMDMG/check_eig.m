function [full_A, eig_value, sym_err] = check_eig(DegDM,BlockL,SupMat,SubMat,if_plot)
size_L = size(BlockL);

N = size_L(1)*size_L(3); % size of A
full_A = zeros(N,N);
for i = 1:N
    e_i = zeros(N,1); e_i(i) = 1;
    full_A(:,i) = MLx(e_i,DegDM,BlockL,SupMat,SubMat);
end

eig_value = eig(full_A);
sym_err = max(max((full_A-full_A')));

if if_plot==1
    figure(5)
    plot(eig_value,zeros(size(eig_value)),'o')
    title('Eigenvalues spectrum of matrix')
    
    figure(6)
    spy(full_A)
    title(['spy of matrix, max(abs(A^T-A)) = ' num2str(sym_err)])
end


end