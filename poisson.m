function [u,x,y] = poisson(f,g,bx0,bxf,by0,byf,D,Mx,My,tol,MaxIter)
% solve u_xx + u_yy + g(x,y)u = f(x,Y)
% over the region D=[x0,xf,y0,yf] = {(x,y)|x0<=x<=xf , y0<= y <= yf}
% u(x0,y)=bx0(y),u(xf,y)=bxf(y)
% u(x,y0)=by0(x),u(x,yf)=byf(x)
% Mx = # of subintervals along x axis
% My = # of subintervals along y axis
% tol     : error tolerance



% MaxIter : the maximum # of iterations
x0 = D(1); xf = D(2); y0 = D(3); yf = D(4);%no need to define?
dx = (xf - x0)/Mx; x = x0 + [0:Mx]*dx;
dy = (yf - y0)/My; y = y0 + [0:My]*dy;
Mx1 = Mx + 1; My1 = My + 1;

%---------------------------------------------------
u0 = zeros(My1, Mx1); % old or the initial guess of soulution 
u  = zeros(My1, Mx1); % new solution
F  = zeros(My,  Mx);
G  = zeros(My,  Mx);
%---------------------------------------------------

%Boundary conditions
for m=1:My1
  u(m,[1 Mx1])=[bx0(y(m)) bxf(y(m))]; 
end %left/right side

for n=1:Mx1
  u([1 My1],n)=[by0(x(n)) byf(x(n))]; 
end %bottom/top

%initialize as the average of boundary values
%sum_of_bv = sum(sum([u(2:My,[1 Mx1]) u([1 My1],2:Mx)']));
%u0(2:My,2:Mx) = sum_of_bv/(2*(Mx + My - 2));

%---------------------------------------------------
% set the initial solution
u0(2:My,2:Mx) = 0.0; 
%---------------------------------------------------

for i = 1:My
  for j = 1:Mx
    F(i,j) = f(x(j),y(i));
    G(i,j) = g(x(j),y(i));
  end
end

dx2  = dx*dx;
dy2  = dy*dy;
dxy2 = 2*(dx2+dy2);
rx   = dx2/dxy2; ry = dy2/dxy2; rxy = rx*dy2;

% used to record the convergence history
itrs  = [];
errs = [];

% update the new solution with u0 using Eq.(9.1.5a)
for itr = 1:MaxIter
  for j = 2:Mx
    for i = 2:My
      u(i,j) = ry *(u0(i,j+1)+u0(i,j-1))+rx*(u0(i+1,j)+u0(i-1,j)) + ...
               rxy*(G(i,j)  *u0(i,j)-F(i,j)); %Eq.(9.1.5a)
    end
  end
  
  err = max(max(abs(u-u0))); % the max value of err
  itrs = [itrs, itr];
  errs = [errs, err];
  
  if itr > 1 && err < tol
    break;
  end
  
  u0 = u;
  
  if mod(itr, 100) == 0
    figure(1); clf;
    subplot(1,2,1);
    surfc(x,y,u); colorbar;
    axis([0 4 0 4]);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(['solution, iteration = ', num2str(itr, '%d'), ', error = ', num2str(err, '%g')]);
    
    subplot(1,2,2);
    semilogy(itrs, errs);
    grid on;
    xlabel('iteration');
    ylabel('error');
    title('convergence history');
    
    drawnow;
  end
end