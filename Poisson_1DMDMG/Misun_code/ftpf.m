%    FTP data to mcs

function ftpf_put(string_in)

mcs = ftp('ftp.mcs.anl.gov');
cd(mcs,'pub/incoming');
mput(mcs,string_in);
