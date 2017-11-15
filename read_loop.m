function [loop,info] = read_loop(cname,lname)
% read .loop file to update parameters in each cases

fpath     = ['cases/' cname];
lpath     = ['cases/' cname '/' lname '.loop'];

fid = fopen(lpath);

tline = fgets(fid);
Nloop = sscanf(tline,'%d');

tline = fgets(fid);
Nvar  = sscanf(tline,'%d');

tline = fgets(fid);
lmode = sscanf(tline,'%d');

info.Nloop = Nloop; % size of loops
info.mode  = lmode; % changing mode, 0=param_i, 1=eval
info.Nvar  = Nvar;  % number of changing variables

switch lmode
    case 0
        tline = fgets(fid);
        for il = 1:Nloop
            tline = fgets(fid);
            tlines= strsplit(tline,',');
            ind   = sscanf(tlines{1},'%d',[1,1]);
            str   = sscanf(tlines{2},' %s ');
            for iv = 1:Nvar
                vname(iv) = sscanf(tlines{2*iv+1},'%d',[1,1]);
                vvalue(iv)= sscanf(tlines{2*iv+2},'%f',[1,1]);
            end
            loop(ind).name = str;
            loop(ind).param_i   = vname;
            loop(ind).var_value = vvalue;
        end
    case 1
        tline = fgets(fid);
        for il = 1:Nloop
            tline = fgets(fid);
            tlines= strsplit(tline,',');
            ind   = sscanf(tlines{1},'%d',[1,1]);
            str   = sscanf(tlines{2},' %s ');
            for iv = 1:Nvar
                vname{iv} = sscanf(tlines{2*iv+1},' %s ');
                vvalue(iv)= sscanf(tlines{2*iv+2},'%f',[1,1]);
            end
            loop(ind).name = str;
            loop(ind).var_name = vname;
            loop(ind).var_value= vvalue;
        end
    case 2
        tline = fgets(fid);

        tline = fgets(fid);
        tlines= strsplit(tline,',');
        ind   = sscanf(tlines{1},'%d',[1,1]);
        str   = sscanf(tlines{2},' %s ');
        for il = 1:Nloop
            for iv = 1:Nvar
                vname{iv} = sscanf(tlines{2*iv+1},' %s ');
                veval{iv}= sscanf(tlines{2*iv+2},' %s ');
            end
            loop(il).name = str;
            loop(il).var_name = vname;
            loop(il).var_eval = veval;
        end
end

fclose(fid);
end
