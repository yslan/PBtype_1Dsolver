
clear all

% ========================= START of PARAMETER ========================= %
% Info of Domain
    NN = 10;
    TotNumDM = 32;
    Nc = 4;%floor(NN/2); % coarse grid size

    xmin = 1; 
    xmax = 2.0;
% Ploting control
solve_plot = 1; % all-in-one plot.

fun_case = 1; % chosing test case
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

%   transform f_fun/a_fun to f_fun_vec/a_fun_vec
    Xendpt = xmin + (xmax-xmin)/TotNumDM*(0:TotNumDM);
    DegDM = NN * ones([1 TotNumDM]);
    [Xdomain, ~] = MDGrid(DegDM,Xendpt);
    Xprolong = StackVec(DegDM,Xdomain);
    
    f_fun_vec = f_fun(Xprolong);
    a_fun_vec = a_fun(Xprolong);
    
    xm = Xprolong(1);   % domain left-end point 
    xp = Xprolong(end); % domain right-end point
    u_endpt = [u_ex(xm),u_ex(xp)];
    dudx_endpt = [dudx(xm),dudx(xp)];

[vlong_usevc,iter_usevc,label_usevc,...
        x_vc,label_vc] =...
    Poisson_1DMDMG(NN,TotNumDM,Nc,xmin,xmax,...
    a_fun_vec,f_fun_vec,u_endpt,dudx_endpt);



% PART 2: solve_plot
if solve_plot==1
    figure(2)
    
    subplot(2,1,1)       % add second plot in 2 x 2 grid
    p_ex  = plot(Xprolong,u_ex(Xprolong),'b-','MarkerSize',14);hold on
    p_sol = plot(Xprolong,vlong_usevc,'r.','MarkerSize',14);
    title({['vcycle, #iter =' num2str(iter_usevc) ',   NN =' num2str(NN)];...
           ['Solver: ' label_usevc]})
    legend([p_ex,p_sol],'exact','solved','Location','southeast')

    subplot(2,1,2)       % add fourth plot in 2 x 2 grid
    err_usevc = abs(u_ex(Xprolong)- vlong_usevc);
    semilogy(Xprolong,err_usevc,'.','MarkerSize',14)
    title({['error: using vcycle, #iter =' num2str(iter_usevc) ',   NN =' num2str(NN)];...
           ['Solver: ' label_usevc ', max error =' num2str(max(err_usevc))]})
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
        end
end