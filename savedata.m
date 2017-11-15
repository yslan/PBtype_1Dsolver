function savedata(cname,filename,index,X,C1,C2,Phi)
% save data into <filename>_<index>.csv
% if index = 0, no index on file name

Nx = length(X);

sind = sprintf('%03d', index);
if index ~=0
    opath = ['cases/' cname '/' filename '_' sind '.csv'];
else
    opath = ['cases/' cname '/' filename '.csv'];
end
fid = fopen(opath,'w');

fprintf(fid,'X,C1,C2,Phi\n');
for ii = 1:Nx
    fprintf(fid,'%4.6f,%4.6f,%4.6f,%4.6f \n',X(ii),C1(ii),C2(ii),Phi(ii));
end

fclose(fid);

end
