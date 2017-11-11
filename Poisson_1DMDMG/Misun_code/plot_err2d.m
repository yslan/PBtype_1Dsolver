function[ip] = plot_err(u,ex,X,Y,Vx,Vy,Lx,Ly,str,v,ihold,ispec); 

ee=ex-u; sp=Vx'*ee*Vy;  %% Spectrum of error

norm(ee)



hold off; if ihold==1; hold on; end; 

if ispec==0; 
       [nx,ny]=size(ee); Px=speye(nx); Py=speye(ny);
       Px=[zeros(1,nx); Px ; zeros(1,nx)];
       Py=[zeros(1,ny); Py ; zeros(1,ny)];
       ee=Px*ee*Py'; % Prolongate error to boundary, for plotting.
       mesh(X,Y,ee);
       strg=sprintf('%s','Error after ',int2str(v),' V-cycles');
       xlabel('x'); ylabel('Solution / Error'); title(strg); axis equal
else;
       hold off; mesh(sp); hold on;
       mesh(.01+0*sp); view([1 1 1])
       strg=sprintf('%s','Spectrum after ',int2str(v),' V-cycles');
       xlabel('x'); ylabel('Spectrum / Error'); title(strg); axis square
end;

ip=1;
if ihold==1; pause; end;

