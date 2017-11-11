function read_rea(cname,rname)
% read rea file
global param

reapath   = ['cases/' cname '/' rname '.rea'];

fid = fopen(reapath,'r');

tline = fgets(fid); tline = fgets(fid);
Nparam = sscanf(tline,'%d');

param = zeros(Nparam,1);
for i=1:Nparam
   tline = fgets(fid);
   param(i) = sscanf(tline,'%f',[1,1]);
end

fclose(fid);
end
