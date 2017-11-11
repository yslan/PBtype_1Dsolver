function [vlong_usevc,iter_usevc,label_usevc,...
        x_vc, label_vc] =...
    Poisson_1DMDMG(NN,TotNumDM,Nc,xmin,xmax,...
    a_fun_vec,f_fun_vec,u_endpt,dudx_endpt)

% Choose linear solver, solving Ax=b
%   inside the v-cycle
    sType_vc = 1;  % CG
    sType_fin = 1; % CG
% needed in any Iterative method
    tol = 1e-6;
    maxit = 6000;
    xc_inig = 0.5*ones((Nc+1)*TotNumDM,1); % initial guess of iter. method
    x_inig = 0.5*ones((NN+1)*TotNumDM,1);  % initial guess of iter. method
% Choose Algorithm in the smoothing part of V-clcle
    algo_vc = 1; % Non-overlapped Jacobi Schwartz


% Construct linear system Lu = Fprolong, on grid Xprolong
%   on fine grid
[Xprolong,DegDM,BlockL,SupMat,SubMat,Fprolong] =...
            get_L_Fvec(xmin,xmax,NN,TotNumDM,u_endpt,dudx_endpt,a_fun_vec,f_fun_vec);
inv_BlockL = cal_inv_BlockL(TotNumDM,BlockL); % pre-compute inverse
%   on coarse grid
%   Construct projection J, from coarse grid into fine grid
    Xendpt   = xmin + (xmax-xmin)/TotNumDM*(0:TotNumDM);
    DegDMc=Nc * ones([1 TotNumDM]);
    [Xcdomain, ~] = MDGrid(DegDMc,Xendpt);
    Xcprolong = StackVec(DegDMc,Xcdomain);

    addpath([pwd '/Misun_code']) % use Misun's code
    Ie = speye(TotNumDM);
    Jh = interp_mat(Xprolong(1:NN+1),Xcprolong(1:Nc+1));
    JL = kron(Ie,Jh);
%   transform f_fun to F_fun_vec
    fc_fun_vec = interp1q(Xcprolong,f_fun_vec,Xcprolong); % useless info.
    ac_fun_vec = interp1q(Xcprolong,a_fun_vec,Xcprolong);
%     fc_fun_vec = JL'*f_fun_vec;
%     ac_fun_vec = JL'*a_fun_vec;
[Xcprolong,DegDMc,BlockLc,SupMatc,SubMatc,Fcprolong] =...
            get_L_Fvec(xmin,xmax,Nc,TotNumDM,u_endpt,dudx_endpt,ac_fun_vec,fc_fun_vec);

% Construct pre-condition matrix of Jacobi method
tmp = zeros((NN+1),TotNumDM); % pre-condition M, Jacobi
for i =1:TotNumDM
    tmp(:,i) = diag(BlockL(:,:,i));
end
M = reshape(tmp,(NN+1)*TotNumDM,1); M = diag(1./M);


% V-cycle
%   settings
N_vcycle = 10;
m_smooth = 4;
sigma = 2/3; % sigma =  .3; % sigma = input('Input sigma:');
%   initializing
r = Fprolong;                               % RHS, initial residue
b = r; 
x_vc_ini = x_inig;     % initial guess
for vcycle = 1:N_vcycle;
	x_vc = x_vc_ini; 
	r = b-MLx(x_vc,DegDM,BlockL,SupMat,SubMat); % residue
    
    % smoothing
    switch algo_vc % Choose Smoothing Algo.
        case 0 % Jacobi smoothing
            for k=1:m_smooth;
                x_vc = x_vc + sigma*M*r;
                r = b - MLx(x_vc,DegDM,BlockL,SupMat,SubMat);
            end
        case 1 % Additive Jacobi Schwartz
            for k = 1:m_smooth
                z = 0*x_vc;
                ind_0 =1; ind_i = NN+1;
%                 [z, iter, label] = ...
%                         choose_solver(TotNumDM,4,tol,maxit,0*z,r,BlockL,DegDM,0*SupMat,0*SubMat);
                for ind_DM = 1:TotNumDM
                    z_sol = inv_BlockL(:,:,ind_DM)*r(ind_0:ind_i);
                    z(ind_0:ind_i) = z(ind_0:ind_i) + z_sol;
                    ind_0 = ind_0 + NN+1; ind_i = ind_i + NN+1;
                end
                x_vc = x_vc + sigma*z;
                r = b - MLx(x_vc,DegDM,BlockL,SupMat,SubMat);
            end
        case 2 % Overlapped Jacobi Schwartz 
            for k = 1:m_smooth
                z = 0*x_vc;
                ind_0 =1; ind_i = NN+1;
                
                % Start of the Loop
                %   The Left-Top Block
                    ind_DM = 1;
                    B_inv = inv_block_mat(inv_BlockL(:,:,ind_DM),...
                            SupMat(:,1,ind_DM),BlockL(1,1,ind_DM+1),1);
                    z_sol = B_inv*r(ind_0:ind_i+1);
                    z(ind_0:ind_i+1) = z(ind_0:ind_i+1) + z_sol;
                    ind_0 = ind_0 + NN+1; ind_i = ind_i + NN+1;
                %   Intermediate Blocks
                for ind_DM = 2:TotNumDM-1
                    B_inv = inv_block_mat(inv_BlockL(:,:,ind_DM),...
                            SupMat(:,1,ind_DM),BlockL(1,1,ind_DM+1),1);
                    B_inv_2 = inv_block_mat(B_inv,...
                            [0;SubMat(:,end,ind_DM-1)],BlockL(end,end,ind_DM-1),2);
                    z_sol = B_inv_2*r(ind_0-1:ind_i+1);
                    z(ind_0-1:ind_i+1) = z(ind_0-1:ind_i+1) + z_sol;
                    ind_0 = ind_0 + NN+1; ind_i = ind_i + NN+1;
                end
                %   The Right-Bottom Block
                    ind_DM = TotNumDM;
                    B_inv = inv_block_mat(inv_BlockL(:,:,ind_DM),...
                            SubMat(:,end,ind_DM-1),inv_BlockL(end,end,ind_DM-1),2);
                    z_sol = B_inv*r(ind_0-1:ind_i);
                    z(ind_0-1:ind_i) = z(ind_0-1:ind_i) + z_sol;
                    ind_0 = ind_0 + NN+1; ind_i = ind_i + NN+1;
                % End of the Loop
                x_vc = x_vc + sigma*z;
                r = b - MLx(x_vc,DegDM,BlockL,SupMat,SubMat);
            end
    end % end of choosing smoothing Algo.
    
    % Coarse-grid restriction  x <--- x + e, where e is approximated on coarse grid.
    rc = JL'*r;    % Restrict residual to coarse grid
    %   Solve coarse-grid problem:  Ac*ec = rc;
    [ec,~,label_vc] =...
        choose_solver(TotNumDM,sType_vc,tol,maxit,xc_inig,rc,BlockLc,DegDMc,SupMatc,SubMatc);
	
    ef = JL*ec;    % Interpolate coarse-grid error to fine grid
	x_vc = x_vc + ef;
    
    % treat result as initial in the next loop
	x_vc_ini = x_vc;
end % end of Vcycle


% Solve results
    tic
    [vlong_usevc, iter_usevc, label_usevc] =...
        choose_solver(TotNumDM,sType_fin,tol,maxit,x_vc,Fprolong,BlockL,DegDM,SupMat,SubMat);
    t_usevc = toc;
end






