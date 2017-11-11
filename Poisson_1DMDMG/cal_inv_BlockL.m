function inv_BlockL = cal_inv_blockl(TotNumDM,BlockL)
if TotNumDM==1
    inv_BlockL = BlockL^(-1);
else
    inv_BlockL = zeros(size(BlockL));
    for ind_DM = 1:TotNumDM
        inv_BlockL(:,:,ind_DM) = BlockL(:,:,ind_DM)^(-1);
    end
end



end