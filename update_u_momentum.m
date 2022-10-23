function [u,aP,aE,aW,aN,aS,Su] = update_u_momentum(Lx,Ly,u,ue,uw,v,vn,...
    vs,M,N,mu,rho,p,alpha_u)

Fe = rho * ue;
Fw = rho * uw;
Fn = rho * vn;
Fs = rho * vs;
dx = Lx  / M;
dy = Ly  / N;
Dx = mu  / dx;
Dy = mu  / dy;


x = zeros(M, 1 );
y = zeros(1 , N);

x(1:M,1) = dx/2 + [0:M-1] * dx;
y(1,1:N) = dy/2 + [0:N-1] * dy;

aW = zeros(M,N);
% aW(1:M,1) = 0 ;
% aW(1:M,2:N) = Dx + Fx(1:M,2:N);
for j = 1:N
    for i = 1:M
            aW(i,j) = Dx + max(Fw(i,j),0);
    end
end

aE = zeros(M,N);
% aE(1:M,1:N-1) = Dx;
% aE(2:M, N)  = 0 ;
for i = 1:M
    for j = 1:N
            aE(i,j) = Dx + max(0,-Fe(i,j));
    end
end

aN = zeros(M,N) ;
% aN(1,1:N) = 0;
% aN(2:M,1:N) = Dy -Fy(2:M,1:N);
for i = 1:M
    for j = 1:N
            aN(i,j) = Dx + max(Fs(i,j),0);
    end
end

aS = zeros(M,N) ;
% aS(1:M-1,1:N) = Dy;
% aS(  M  ,1:N) = 0 ;
for i = 1:M
    for j = 1:N
            aS(i,j) = Dx + max(0,-Fn(i,j));
    end
end

SP = zeros(M,N);
SP(2:M-1,1) = -(2*Dx + max(Fw(i,j),0));
SP(2:M-1,N) = -(2*Dx + max(0,-Fe(i,j)));
SP(1,2:N-1) = -(2*Dy + max(0,-Fn(i,j)));
SP(M,2:N-1) = -(2*Dy + max(Fs(i,j),0));

SP( 1 , 1 ) = SP(2,1)+SP(1,2);
SP( 1 , N ) = SP(2,N)+SP(1,N-1);
SP( M , 1 ) = SP(M-1,1)+SP(M,2);
SP( M , N ) = SP(M-1,N)+SP(M,N-1);
% SP( 1 , 1 ) = -1e10;
% SP( 1 , N ) = -1e10;
% SP( M , 1 ) = -1e10;
% SP( M , N ) = -1e10;

aP = aW + aE + aS + aN -SP + Fe - Fw - Fn + Fs;
% mesh(x,y,aP);

Su  = zeros(M,N);
% Suu = zeros(M,N);
% for j = 2:N-1 
%     Suu(1,j) = (2*Dy + Fy(1,j))*u(1,j);
% end 
% for j = 2:N-1 
%     Suu(M,j) = (2*Dy - Fy(M,j))*u(M,j);
% end 
% for i = 2:M-1 
%     Suu(i,1) = (2*Dx + Fx(i,1))*u(i,1);
% end 
% for i = 2:M-1 
%     Suu(i,N) = (2*Dx - Fx(i,N))*u(i,N);
% end

for j = 2:N-1 
    Su(1,j) = (2*Dy + max(0,-Fn(i,j)))*u(1,j)+ 2 * dy * (p(1,j-1) - p(1,j+1)) / aP(i,j);
end 
for j = 2:N-1 
    Su(M,j) = (Dx + max(Fs(i,j),0))*u(M,j)+ 2 * dy * (p(M,j-1) - p(M,j+1)) / aP(i,j);
end 
for i = 2:M-1 
    Su(i,1) = ((2*Dx + max(Fw(i,j),0)))*u(i,1)+ 2 * dy * (p(i,1) - p(i,2)) / aP(i,j);
end 
for i = 2:M-1 
    Su(i,N) = (2*Dx + max(0,-Fe(i,j)))*u(i,N)+ 2 * dy * (p(i,N-1) - p(i,N)) / aP(i,j);
end
for i = 2:M-1
    for j = 2:N-1
        Su(i,j) = aP(i,j) / alpha_u * (1-alpha_u) * (u(i,j)) + 2 * dy * (p(i,j-1) - p(i,j+1)) / aP(i,j);
    end
end
Su( 1 , 1 ) = (Su(1,2)+Su(2,1));
Su( 1 , N ) = (Su(1,N-1)+Su(2,N));
Su( M , 1 ) = (Su(M-1,1)+Su(M,2));
Su( M , N ) = (Su(M-1,N)+Su(M,N-1));

re_Su = reshape_a(Su,M,N);

%%
%对Su重新编号
aP = aP / alpha_u;
matrix = sparse_coef_auto(aP,aW,aE,aN,aS,M,N);

phi = matrix\re_Su;

re_phi = zeros(M,N);

for j = 1 : N
    for i = 1 : M
        re_phi(i,j) = phi((j-1)*M+i,1);
    end
end
u = re_phi;