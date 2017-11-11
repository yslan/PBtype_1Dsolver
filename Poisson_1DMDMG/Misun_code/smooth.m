

%  Illustration of Damped-Jacobi smoothing effect on error

%  Usage:   N=40; k=3; smooth

z = 0:N; z=z'/N; h = z(2)-z(1); % Points and gridspacing on [0,1]
n = N-1; e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n) / (h*h);
D = diag(diag(A)); sigma=0.66666;    % Smoother

ex= rand(n,1); f = A*ex;   % f = A*exact
ex= [0; ex; 0];  % Extend exact solution to boundaries (for plotting)

u = 0*f;                       % Initial guess
for j=1:k                      % Apply smoother k times
   u = u + sigma*(D\(f-A*u));  % Smoothing step
end;

%DIAGNOSITICS: Use known solution to show error behavior.
ue = [0; u; 0]; ee =ex-ue; plot(z,ue,'bo-',z,ee,'r-',z,ex,'k--',z,0*ex,'k-'); 


strg=sprintf('%s','Error after ',int2str(k),' smoothings');
xlabel('x'); ylabel('Solution / Error'); title(strg); axis square
legend('numeric','error','exact'); norm(ee); pause(0.1)

set (gcf,'color',[1 1 1]) ; figure;
plot(z,ue,'bo-',z,0*ex,'k-'); axis square; axis([0 1 -.5 1.]);
xlabel('x'); ylabel('Current Estimate for u_k'); title(strg); axis square
set (gcf,'color',[1 1 1]) ; figure;
plot(z,ee,'r-' ,z,0*ex,'k-'); axis square; axis([0 1 -.5 1.]);
xlabel('x'); ylabel('Error'); title(strg); axis square
set (gcf,'color',[1 1 1]) ; figure;
plot(z,ex,'k--',z,0*ex,'k-'); axis square; axis([0 1 -.5 1.]);
xlabel('x'); ylabel('Exact Solution'); title(strg); axis square
set (gcf,'color',[1 1 1]) ;

