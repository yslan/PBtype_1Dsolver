function init_loop_folder(cname,rname,lname)
% check if folder rea userf.m exist. If not, create an example

fpath     = ['cases/' cname];
reapath   = ['cases/' cname '/' rname '.rea'];
usrpath   = ['cases/' cname '/userf.m'];
looppath  = ['cases/' cname '/' lname '.loop'];

if exist(fpath, 'dir')==0 % folder not exist
    % create folder, rea, user file, and param
    [status, msg, msgID] = mkdir(fpath);
    st1  = copyfile('cases/loop_test/userf.m',fpath);
    st21 = copyfile('cases/loop_test/test.rea',fpath);
    st22 = movefile([fpath '/test.rea'],reapath);
    st31 = copyfile('cases/loop_test/loop.rea',fpath);
    st32 = movefile([fpath '/loop.rea'],reapath);
    error('There is no such folder, create an example and exit\n');
else
    if exist(reapath,'file') == 0 % rea not exist
    % create example rea
        st21 = copyfile('cases/loop_test/test.rea',fpath);
        st22 = movefile([fpath '/test.rea'],reapath);
        error('There is no .rea, create an example and exit\n');
    elseif exist(usrpath,'file') == 0 % usr.m not exist
    % create example userf.m
        st1 = copyfile('cases/loop_test/userf.m',fpath);
        error('There is no userf.m, create an example and exit\n');
    elseif exist(looppath,'file') == 0 % usr.m not exist
    % create example loop
        st31 = copyfile('cases/loop_test/test_mode0.loop',fpath);
        st32 = movefile([fpath '/test_mode0.loop'],reapath);
        error('There is no .loop, create an example and exit\n');
    end
end
