function init_folder(cname,rname)
% check if folder rea userf.m exist. If not, create an example

fpath     = ['cases/' cname];
reapath   = ['cases/' cname '/' rname '.rea'];
usrpath   = ['cases/' cname '/userf.m'];

if exist(fpath, 'dir')==0 % folder not exist
    % create folder, rea, user file, and param
    [status, msg, msgID] = mkdir(fpath);
    st1 = copyfile('cases/ClassicalPB_1110/userf.m',fpath);
    st2 = copyfile('cases/ClassicalPB_1110/test.rea',fpath);
    st3 = movefile([fpath '/test.rea'],reapath);
    error('There is no such folder, create an example and exit\n');
else
    if exist(reapath,'file') == 0 % rea not exist
    % create example rea
        st2 = copyfile('cases/ClassicalPB_1110/test.rea',fpath);
        st3 = movefile([fpath '/test.rea'],reapath);
        error('There is no .rea, create an example and exit\n');
    elseif exist(usrpath,'file') == 0 % usr.m not exist
    % create example userf.m
        st1 = copyfile('cases/ClassicalPB_1110/userf.m',fpath);
        error('There is no userf.m, create an example and exit\n');
    end
end
