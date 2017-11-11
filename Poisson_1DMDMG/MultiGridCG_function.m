function MultiGridCG_function(parameter,save_info)
% This function is used to test the Multi-Grid(MG) algorithm

    fig_save = save_info.fig_save;
    mat_save = save_info.mat_save;
    out_file_path = save_info.out_file_path;
    save_label = save_info.save_label;
    
    NN = parameter.NN;
    TotNumDM = parameter.TotNumDM;
    Nc = parameter.Nc;
    xmin = parameter.xmin;
    xmax = parameter.xmax;
    fun_case = parameter.fun_case;
    sType_vc = parameter.sType_vc;
    sType_fin = parameter.sType_fin;
    algo_vc = parameter.algo_vc;
    tol = parameter.tol;
    maxit = parameter.maxit;
    N_vcycle = parameter.N_vcycle;
    m_smooth = parameter.m_smooth;
    sigma = parameter.sigma;

% Ploting control
vc_plot = 1; % 1: plot error of v-cycle, 2: animation, on fig. 1
solve_plot = 1; % all-in-one plot.

xc_inig = 0.5*rand((Nc+1)*TotNumDM,1); % initial guess of iter. method
x_inig = 0.5*rand((NN+1)*TotNumDM,1);  % initial guess of iter. method
switch fun_case
    case 0 % bounday layer
        epsilon = 1/10/2;
        u_ex = @(x) -(2*x-3).^(1/epsilon) +1 ; % sin(4*pi*x);
        dudx = @(x) -(2/epsilon) * (2*x-3).^(1/epsilon-1); %4*pi*cos(4*pi*x);
        a_fun = @(x) epsilon*x; 
        f_fun = @(x) 2 *( (2*x-3).^(1/epsilon-1) ...
                          + 2 *(1/epsilon-1) * (x .* (2*x-3).^(1/epsilon-2)));
    case 1 % smooth function
        u_ex = @(x) 1+cos(pi*x);
        dudx = @(x) -pi*sin(pi*x);
        a_fun = @(x) ones(size(x));
        f_fun = @(x) +pi^2*cos(pi*x);
    case 2 % discontinuity function
        domain_plot = 1; % plot or not 
        if mod(TotNumDM,2) == 1
            error('ERR: In this test case, TotNumDM must be even. test case: 2'); 
        end

        m = 1; %int
        n = 2; %int
        x_mid = (xmax+xmin)/2;
        step_fun_r = @(x) heaviside(x-x_mid);
        step_fun_l = @(x) 1-heaviside(x-x_mid);
        
        u_1     = @(x) 1+sin(m*pi*(x-x_mid));
        dudx_1  = @(x) m*pi*cos(m*pi*(x-x_mid));
        a_fun_1 = @(x) ones(size(x));
        f_fun_1 = @(x) +(m*pi)^2*sin(m*pi*(x-x_mid));
        u_2     = @(x) 1+sin(n*pi*(x-x_mid));
        dudx_2  = @(x) n*pi*cos(n*pi*(x-x_mid));
        a_fun_2 = @(x) m/n*ones(size(x));
        f_fun_2 = @(x) +(n*m)*pi^2*sin(n*pi*(x-x_mid));
        
        u_ex  = @(x) u_1(x)    .*step_fun_l(x) + u_2(x)    .*step_fun_r(x);
        dudx  = @(x) dudx_1(x) .*step_fun_l(x) + dudx_2(x) .*step_fun_r(x);
        a_fun = @(x) a_fun_1(x).*step_fun_l(x) + a_fun_2(x).*step_fun_r(x); 
        f_fun = @(x) f_fun_1(x).*step_fun_l(x) + f_fun_2(x).*step_fun_r(x);
end
% ========================== END of PARAMETER ========================== %




% Construct linear system Lu = Fprolong, on grid Xprolong
%   on fine grid
[Xprolong,DegDM,BlockL,SupMat,SubMat,Fprolong] =...
            get_L(xmin,xmax,NN,TotNumDM,u_ex,dudx,a_fun,f_fun);
inv_BlockL = cal_inv_BlockL(TotNumDM,BlockL); % pre-compute inverse
%   on coarse grid
[Xcprolong,DegDMc,BlockLc,SupMatc,SubMatc,Fcprolong] =...
            get_L(xmin,xmax,Nc,TotNumDM,u_ex,dudx,a_fun,f_fun);

% Construct projection J, from coarse grid into fine grid
addpath([pwd '/Misun_code']) % use Misun's code
Ie = speye(TotNumDM);
Jh = interp_mat(Xprolong(1:NN+1),Xcprolong(1:Nc+1));
JL = kron(Ie,Jh);

Jh_f2c = interp_mat(Xcprolong(1:Nc+1),Xprolong(1:NN+1)); % testing
JL_f2c = kron(Ie,Jh_f2c);

% Construct pre-condition matrix of Jacobi method
tmp = zeros((NN+1),TotNumDM); % pre-condition M, Jacobi
for i =1:TotNumDM
    tmp(:,i) = diag(BlockL(:,:,i));
end
M = reshape(tmp,(NN+1)*TotNumDM,1); M = diag(1./M);




% V-cycle
%   storage
Xvc_smooth_store = zeros(length(Xprolong),N_vcycle);
Xvc_recorrect_store = zeros(length(Xprolong),N_vcycle);
%   initializing
r = Fprolong;                               % RHS, initial residue
b = r; 
x_vc_ini = 0.5*rand((NN+1)*TotNumDM,1);     % initial guess
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
    Xvc_smooth_store(:,vcycle) = x_vc; % store results

    
    % Coarse-grid restriction  x <--- x + e, where e is approximated on coarse grid.
    rc = JL'*r;    % Restrict residual to coarse grid
    %   Solve coarse-grid problem:  Ac*ec = rc;
    if sType_vc == 0
        [ec, iter_vc, label_vc] =...
            choose_solver(TotNumDM,sType_vc,tol,maxit,xc_inig,rc,BlockLc);
    else 
        [ec, iter_vc, label_vc] =...
            choose_solver(TotNumDM,sType_vc,tol,maxit,xc_inig,rc,BlockLc,DegDMc,SupMatc,SubMatc);
    end
	ef = JL*ec;    % Interpolate coarse-grid error to fine grid
	x_vc = x_vc + ef;
    Xvc_recorrect_store(:,vcycle) = x_vc; % store results
    
    % treat result as initial in the next loop
	x_vc_ini = x_vc;
end; % end of Vcycle




% Solve results
if sType_fin == 0
    tic
    [vlong_usevc, iter_usevc, label_usevc] =...
        choose_solver(TotNumDM,sType_fin,tol,maxit,x_vc,Fprolong,BlockL);
    t_usevc = toc
    tic
    [vlong_direct, iter_direct, label_direct] =...
        choose_solver(TotNumDM,sType_fin,tol,maxit,x_inig,Fprolong,BlockL);
    t_direct = toc
else 
    tic
    [vlong_usevc, iter_usevc, label_usevc] =...
        choose_solver(TotNumDM,sType_fin,tol,maxit,x_vc,Fprolong,BlockL,DegDM,SupMat,SubMat);
    t_usevc = toc
    tic
    [vlong_direct, iter_direct, label_direct] =...
        choose_solver(TotNumDM,sType_fin,tol,maxit,x_inig,Fprolong,BlockL,DegDM,SupMat,SubMat);
    t_direct = toc
end




% Plot result
% PART 1: vc_plot
switch vc_plot
    case 1
        figure(1)
        max_eL_s = zeros(N_vcycle,1);
        max_eL_r = zeros(N_vcycle,1);
        for k = 1:N_vcycle
            max_eL_s(k) = max(abs(Xvc_smooth_store(:,k) - u_ex(Xprolong)));
            max_eL_r(k) = max(abs(Xvc_recorrect_store(:,k) - u_ex(Xprolong)));
        end
            max_ps = semilogy(1:N_vcycle,max_eL_s,'ro'); hold on;
            max_pr = semilogy(1:N_vcycle,max_eL_r,'b-');  
            title({['Error in vcycle, max abs error, TotNumDM =' num2str(TotNumDM)];...
                   ['vcycle solver: ' label_vc ', NN =' num2str(NN)...
                        ', vc final error= ' num2str(max_eL_r(N_vcycle))]})
            legend([max_ps max_pr],'smooth part','restriction part'); 
            hold off;
        if fig_save==1
            fff = gcf;
            fin_save_label = [out_file_path save_label '_fig_1.png'];
            saveas(fff,fin_save_label)
            fin_save_label = [out_file_path save_label '_fig_1.fig'];
            saveas(fff,fin_save_label)
            close('all')
        end
    case 2 % animation
        figure(2)
        for k = 1:N_vcycle
            eL_s = abs(Xvc_smooth_store(:,k) - u_ex(Xprolong));
            eL_r = abs(Xvc_recorrect_store(:,k) - u_ex(Xprolong));
            
            ps = semilogy(Xprolong,eL_s,'ro'); hold on;
            pr = semilogy(Xprolong,eL_r,'b-');  
            title({['Error in vcycle, animation, #cycle = ' num2str(k) ', TotNumDM =' num2str(TotNumDM)];...
                   ['vcycle solver: ' label_vc ', NN =' num2str(NN)...
                        ', vc final error= ' num2str(max_eL_r(N_vcycle))]})
            legend([ps pr],'smooth part','restriction part'); 
            hold off;
            pause
        end
end
% end of PART 1: vc_plot

% PART 2: solve_plot
if solve_plot==1
    figure(21+sType_fin*10)
    suptitle(['Comparing results with vcycle or not, TotNumDM =' num2str(TotNumDM)])
    
    positionVector1 = [0.1, 0.55, 0.4, 0.3];
    subplot('Position',positionVector1)      % add first plot in 2 x 2 grid
    p_ex  = plot(Xprolong,u_ex(Xprolong),'b-','MarkerSize',14);hold on
    p_sol = plot(Xprolong,vlong_direct,'r.','MarkerSize',14);
    title({['directly #iter =' num2str(iter_direct) ',   NN =' num2str(NN)];...
           ['Solver: ' label_direct]})
    legend([p_ex,p_sol],'exact','solved','Location','southeast')

    positionVector2 = [0.55, 0.55, 0.4, 0.3];
    subplot('Position',positionVector2)       % add second plot in 2 x 2 grid
    p_ex  = plot(Xprolong,u_ex(Xprolong),'b-','MarkerSize',14);hold on
    p_sol = plot(Xprolong,vlong_usevc,'r.','MarkerSize',14);
    title({['vcycle, #iter =' num2str(iter_usevc) ',   NN =' num2str(NN)];...
           ['Solver: ' label_usevc]})
    legend([p_ex,p_sol],'exact','solved','Location','southeast')

    positionVector3 = [0.1, 0.1, 0.4, 0.3];
    subplot('Position',positionVector3)       % add third plot in 2 x 2 grid
    err_direct = abs(u_ex(Xprolong)- vlong_direct);
    semilogy(Xprolong,err_direct,'.','MarkerSize',14)
    title({['error: solving directly, #iter =' num2str(iter_direct) ',   NN =' num2str(NN)];...
           ['Solver: ' label_direct ', max error =' num2str(max(err_direct))]})
    
    positionVector4 = [0.55, 0.1, 0.4, 0.3];
    subplot('Position',positionVector4)       % add fourth plot in 2 x 2 grid
    err_usevc = abs(u_ex(Xprolong)- vlong_usevc);
    semilogy(Xprolong,err_usevc,'.','MarkerSize',14)
    title({['error: using vcycle, #iter =' num2str(iter_usevc) ',   NN =' num2str(NN)];...
           ['Solver: ' label_usevc ', max error =' num2str(max(err_usevc))]})
    
   screen_vec = get(0, 'Screensize');
   screen_vec(1) = 10; screen_vec(2) = 50;
   screen_vec(3) = screen_vec(3) - 20; screen_vec(4) = screen_vec(4) -180;
   set(gcf, 'Position', screen_vec);
   
   if fig_save==1
    	fff = gcf;
    	fin_save_label = [out_file_path save_label '_fig_' num2str(21+sType_fin*10) '.png'];
        saveas(fff,fin_save_label)
    	fin_save_label = [out_file_path save_label '_fig_' num2str(21+sType_fin*10) '.fig'];
        saveas(fff,fin_save_label)
        close('all')
    end
end

% PART 3: plot about exact function
switch fun_case
    case 2
        if domain_plot==1
            figure(7)
            x_mid = (xmax+xmin)/2;
            BD_cond_in = a_fun_1(x_mid)*dudx_1(x_mid) - a_fun_2(x_mid)*dudx_2(x_mid);
            plot(Xprolong,dudx(Xprolong))
            title({['dudx'];...
                   ['a_1*dudx_1 - a_2*dudx_2 = ' num2str(BD_cond_in)]})
            if fig_save==1
                fff = gcf;
                fin_save_label = [out_file_path save_label '_fig_7.png'];
                saveas(fff,fin_save_label)
                fin_save_label = [out_file_path save_label '_fig_7.fig'];
                saveas(fff,fin_save_label)
                close('all')
            end
        end
end

% PART 4: save all data into '.dat'
if mat_save==1
    fin_save_label = [out_file_path save_label '.mat'];
    save(fin_save_label)
end


end