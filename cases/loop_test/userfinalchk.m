function [V0,C_D] = userfinalchk(Q,loop,linfo)
% final check after loop, do whatevver you want
global param

global TotNumDM NN xmin xmax
global k_B T e zval epsilon
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff


Nloop = linfo.Nloop;
ND    = Nloop-2; % central difference

V0      = zeros(ND,1);
C_D     = zeros(ND,1);
ex_C_D  = zeros(ND,1);
err_C_D = zeros(ND,1);

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

z = abs(zval(1));
ex_C_D = 2 * sqrt(epsilon*c0*e^2*z^2/(2*k_B*T)) * cosh(e*z*V0/(2*k_B*T));

err_C_D = ex_C_D - C_D;

disp('V0, C_D, ex_C_D,err_C_D')
disp([V0',C_D',ex_C_D',err_C_D'])


figure(10)
p    = plot(V0,C_D,'ro');hold on
p_ex = plot(V0,ex_C_D,'k-');
xlabel('V_0')
legend([p,p_ex],'C_D','C_D exact')

