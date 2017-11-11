
position_down = [0.1 0.1 0.8 0.3];
position_up = [0.1 0.5 0.8 0.3];

% TotNumDM = 128, vcycle = 10, 
xx = [1 2 3 4]; % m_smooth 

c2 = [2379 2355 246 10];
c4 = [2471 2441 216 6];

direct2 = [2354 2354 2354 2354];
direct4 = [2444 2444 2444 2444];

figure(5)
p1 = semilogy(xx,c2,'ro-'); hold on
p2 = semilogy(xx,direct2,'ro--');

p3 = semilogy(xx,c4,'ko-');
p4 = semilogy(xx,direct4,'ko--');
title('m smooth v.s. # of iter, TotNumDM = 128, vcycle = 10')
legend([p1 p2 p3 p4],...
    'use vcycle, c sigma = 2','direct, c sigma = 2',...
    'use vcycle, c sigma = 4','direct, c sigma = 4',...
    'Location','southwest')
xlabel('m smooth')

% m_smooth = 3, vcycle = 10,
figure(1)
xx = [2 3 4]; % c_sigma
direct16 = [262 269 277];
uvc16 = [92 96 98];

suptitle({'c sigma v.s. # of iter, TotNumDM = 16';...
    'm smooth = 3, vcycle = 10'})
subplot('position',position_down)
p1 = plot(xx,direct16,'bo-'); 
subplot('position',position_up)
xlabel('c sigma')
p2 = plot(xx,uvc16,'ro-');
legend([p1 p2],'direct','use vcycle','Location','southeast')
xlabel('c sigma')

figure(2)
direct32 = [558 571 583];
uvc32 = [128 133 134];
suptitle({'c sigma v.s. # of iter, TotNumDM = 32';...
    'm smooth = 3, vcycle = 10'})
subplot('position',position_down)
p1 = plot(xx,direct32,'bo-');
subplot('position',position_up)
xlabel('c sigma')
p2 = plot(xx,uvc32,'ro-');
legend([p1 p2],'direct','use vcycle','Location','southeast')
xlabel('c sigma')


figure(3)
direct64 = [1155 1177 1202];
uvc64 = [189 156 156];
suptitle({'c sigma v.s. # of iter, TotNumDM = 64';...
    'm smooth = 3, vcycle = 10'})
subplot('position',position_down)
p1 = plot(xx,direct64,'bo-');
subplot('position',position_up)
xlabel('c sigma')
p2 = plot(xx,uvc64,'ro-');
legend([p1 p2],'direct','use vcycle','Location','southeast')
xlabel('c sigma')


figure(4)
direct128 = [2352 2400 2444];
uvc128 = [246 212 216];
suptitle({'c sigma v.s. # of iter, TotNumDM = 128';...
    'm smooth = 3, vcycle = 10'})
subplot('position',position_down)
p1 = plot(xx,direct128,'bo-');
subplot('position',position_up)
xlabel('c sigma')
p2 = plot(xx,uvc128,'ro-');
legend([p1 p2],'direct','use vcycle','Location','southeast')
xlabel('c sigma')



% c_sigma = 4, m_smooth = 3,
figure(6)
xx = [16,32,64,128]; % TotNumDM
cg_direct = [378 881 1995 4347];
cg_uv = [35 50 66 104];

gmres_direct = [176 352 704 1408];
gmres_uv = [38 58 83 139];

suptitle({'TotNumDM v.s. # of iter. func = boundary layer';...
    'c sigma = 4, m smooth = 3'})
subplot('position',position_up)
p1 = loglog(xx,cg_uv,'bo-'); hold on
p2 = loglog(xx,gmres_uv,'ko-');
legend([p1 p2],'cg use vcycle','gmres use vcycle','Location','southeast')
xlabel('TotNumDM')

subplot('position',position_down)
p3 = loglog(xx,cg_direct,'bo-'); hold on
p4 = loglog(xx,gmres_direct,'ko-');
legend([p3 p4],'cg direct','gmres direct','Location','southeast')
xlabel('TotNumDM')

% non-overlapping v.s. overlapping
xx = [32 64 128 256 512 1024 2048];
vec_n = [119 174 237 292 399 399 686];
vec_o = [146 193 199 320 338 573 1096];
vec_direct = [2015 4236 8976 18281 36721 73669 147414];

figure(10)
suptitle({'Non-overlapping v.s. overlapping on discontinouous function';
    'NN = 20, Nc = 4, m smooth = 4, vcycle = 10, c sigma = 3'})
subplot('position',position_up)
p1 = semilogx(xx,vec_n,'bo-'); hold on
p2 = semilogx(xx,vec_o,'ro-');
legend([p1 p2],'Non-overlapping in vcycle','Overlapping in vcycle','Location','southeast')
xlabel('TotNumDM')

subplot('position',position_down)
p1 = loglog(xx,vec_n,'bo-'); hold on
p2 = loglog(xx,vec_o,'ro-');
p3 = loglog(xx,vec_direct,'ko-');
legend([p1 p2 p3],'Non-overlapping in vcycle','Overlapping in vcycle','Solving directly','Location','southeast')
xlabel('TotNumDM')
