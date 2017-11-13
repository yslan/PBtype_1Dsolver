function [C1,C2,fval] = solve_c_by_phi(Phi,C1_ini,C2_ini,mode)
% Given Phi, get C1 and C2
%   mode: 
%       fsolve
%       ptws_NT: point wise Newton

global param

global TotNumDM NN xmin xmax
global k_B T e zval
global c0 phi_L Psi_1 Psi_2
global tol_pot tol_c tol_res

global DegDM Xprolong
global MassVec JacVec Diff


Nx = length(Phi);
switch mode
    case 'fsolve'
        option = optimoptions('fsolve','Display','none');
        
        C_ini = [C1_ini;C2_ini];
        [C_out,fval] = fsolve(@(C_in)res_Concen_eq(Phi,C_in),C_ini,option);
        
        C1 = C_out(1:Nx); C2 = C_out(Nx+1:2*Nx);
    
    case 'ptws_NT_couple'
        Jaceps = 1E-8;
        tolNT  = 1E-8;
        maxitNT= 10;
        alphaNT= 0.8; % relaxation, slow down iteration preventing negative

        NTF = @(C)res_Concen_eq(Phi,C);

        ii = 0;
        fval = NTF([C1_ini;C2_ini]);
        resNT = max(abs(fval));
        C1 = C1_ini;
        C2 = C2_ini;

        % NT iter, solve NTF = 0
        while (resNT > tolNT && ii <= maxitNT)
            ii = ii+1;
 
            % approximate Jac mat
            dNTFdC1 = (NTF([C1+Jaceps;C2])-NTF([C1;C2])) / Jaceps;
            dNTFdC2 = (NTF([C1;C2+Jaceps])-NTF([C1;C2])) / Jaceps;
    
            % compute Jac^-1 since 2*2
            detNTF = dNTFdC1(1:Nx).*dNTFdC2(Nx+1:2*Nx) - dNTFdC1(Nx+1:2*Nx).*dNTFdC2(1:Nx);
    
            inv_F  = zeros(Nx,2,2);
            inv_F(:,1,1) =  dNTFdC2(Nx+1:2*Nx)./detNTF;
            inv_F(:,1,2) = -dNTFdC2(1:Nx)./detNTF;
            inv_F(:,2,1) = -dNTFdC1(Nx+1:2*Nx)./detNTF;
            inv_F(:,2,2) =  dNTFdC1(1:Nx)./detNTF;
            
            % main iteration
            s1k = inv_F(:,1,1).*fval(1:Nx) + inv_F(:,1,2).*fval(Nx+1:2*Nx);
            s2k = inv_F(:,2,1).*fval(1:Nx) + inv_F(:,2,2).*fval(Nx+1:2*Nx);

            C1 = C1 - alphaNT * s1k;
            C2 = C2 - alphaNT * s2k;
    
            % compute residual of NT
            fval = NTF([C1;C2]);
            resNT = max(abs(fval));
        end

    case 'ptws_NT_sep'
        Jaceps = 1E-8;
        tolNT  = 1E-8;
        maxitNT= 20;
        alphaNT= 1; % relaxation, slow down iteration preventing negative

        NTF1 = @(C1)res_Concen_eq(Phi,[C1;C2_ini]);
        NTF2 = @(C2)res_Concen_eq(Phi,[C1_ini;C2]);

        ii = 0;
        fval = NTF1(C1_ini);
        fval1= fval(1:Nx);
        resNT = max(abs(fval1));
        C1 = C1_ini;

        % NT iter, solve NTF1 = 0
        while (resNT > tolNT && ii <= maxitNT)
            ii = ii+1;
 
            % approximate Jac mat
            dNTFdC1 = (NTF1(C1+Jaceps)-NTF1(C1)) / Jaceps;
    
            % compute Jac^-1 since 2*2
            inv_F  =  1./dNTFdC1(1:Nx);
            
            % main iteration
            s1k = inv_F.*fval1;
            C1 = C1 - alphaNT * s1k;
    
            % compute residual of NT
            fval = NTF1(C1);
            fval1= fval(1:Nx);
            resNT = max(abs(fval1));
        end


        % NT iter, solve NTF2 = 0
        ii = 0;
        fval = NTF2(C2_ini);
        fval2= fval(Nx+1:2*Nx);
        resNT = max(abs(fval2));
        C2 = C2_ini;
        while (resNT > tolNT && ii <= maxitNT)
            ii = ii+1;
 
            % approximate Jac mat
            dNTFdC2 = (NTF2(C2+Jaceps)-NTF2(C2)) / Jaceps;
    
            % compute Jac^-1 since 2*2
            inv_F =  1./dNTFdC2(Nx+1:2*Nx);
            
            % main iteration
            s2k = inv_F.*fval2
            C2 = C2 - alphaNT * s2k;
    
            % compute residual of NT
            fval = NTF2(C2);
            fval2= fval(Nx+1:2*Nx);
            resNT = max(abs(fval2));
        end
    fval = [fval1;fval2]; 
end
