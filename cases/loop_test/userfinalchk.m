function [V0,C_D] = userfinalchk(Q,loop,linfo)
% final check after loop, do whatevver you want


Nloop = linfo.Nloop;
ND    = Nloop-2; % central difference


% this is bad...
lmode = linfo.mode;
Nvar  = linfo.Nvar;

for il = 1:Nloop % get V_0
switch lmode
    case 0
        for iv = 1:1
            ip = loop(il).param_i(iv);
            param(ip) = loop(il).var_value(iv);
        end

    case 1
        for iv = 1:1
            eval(['V_0(' num2str(il) ') =' num2str(loop(il).var_value) ';']);
        end
    case 2 
        for iv = 1:1
            eval(['varf = ' loop(il).var_eval{iv} ';'])
            eval(['V_0(' num2str(il) ') =' num2str(varf(il)) ';']);
        end
end
end

V0 =  V_0(2:ND+1);
C_D= -(Q(3:ND+2)-Q(1:ND))./(V_0(3:ND+2) - V_0(1:ND));

disp('V0, C_D')
disp([V0',C_D'])
