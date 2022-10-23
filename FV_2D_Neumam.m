function [re_phi,x,y]=FV_2D_Neumam(Lx,Ly,Tx0,Txf,Ty0,Tyf,M,N,k,A,qx0,qy0,qxf,qyf)

dx=Lx/M;
x=dx/2+[0:M-1]*dx;
dy=Ly/N;
y=dy/2+[0:N-1]*dy;

%aW
aW(1:M,1  )   = 0;
aW(1:M,2:N) = k*A/dx;

%aE
aE(1:M,1:N-1) = k*A/dx;
aE(1:M,  N  ) = 0;

%aN
aN(1  ,1:N) = 0;
aN(2:M,1:N) = k*A/dy;

%aS
aS(1:M-1,1:N) = k*A/dy;
aS(  M  ,1:N) = 0;

%SP
SP([1 M],[1 N]) = -2*k*A/dx-2*k*A/dy;
SP(2:M-1,[1 N]) = -2*k*A/dx;
SP([1 M],2:N-1) = -2*k*A/dy;

%aP
aP = zeros(M,N);
aP = aW + aE + aN + aS - SP;

%Su
Su(  1  ,  1  ) = 2*k*A*(Tx0/dx+Ty0/dy)+qx0*A*dx+qy0*A*dy;
Su(  M  ,  1  ) = 2*k*A*(Tx0/dx+Tyf/dy)+qx0*A*dx+qyf*A*dy;
Su(  1  ,  N  ) = 2*k*A*(Txf/dx+Ty0/dy)+qxf*A*dx+qy0*A*dy;
Su(  M  ,  N  ) = 2*k*A*(Txf/dx+Tyf/dy)+qxf*A*dx+qyf*A*dy;
Su(2:M-1,2:N-1) = 0;
Su(2:M-1,1    ) = 2*k*A*Tx0/dx+qx0*A*dx;
Su(2:M-1,N    ) = 2*k*A*Txf/dx+qxf*A*dx;
Su(1    ,2:N-1) = 2*k*A*Ty0/dy+qy0*A*dy;
Su(M    ,2:N-1) = 2*k*A*Tyf/dy+qyf*A*dy;

%%降维
%对Su重新编号
re_Su = zeros(M*N,1);
for j = 1 : N
    for i = 1 : M
        re_Su((j-1)*M+i,1) = Su(i,j) ;
    end
end

%对aW重新编号
re_aW = zeros(M*N,1) ;
for i = 1 : M
    for j = 1 : N
        re_aW((j-1)*M+i,1) = aW(i,j) ;
    end
end
%对aE重新编号
re_aE = zeros(M*N,1);
for i = 1 : M
    for j = 1 : N
        re_aE((j-1)*M+i,1) = aE(i,j) ;
    end
end
%对aN重新编号
re_aN = zeros(M*N,1);
for i = 1 : M
    for j = 1 : N
        re_aN((j-1)*M+i,1) = aN(i,j) ;
    end
end

%对aS重新编号
re_aS = zeros(M*N,1);
for i = 1 : M
    for j = 1 : N
        re_aS((j-1)*M+i,1) = aS(i,j) ;
    end
end

%对aP重新编号
re_aP = zeros(M*N,1);
for i = 1 : M
    for j = 1 : N
        re_aP((j-1)*M+i,1) = aP(i,j) ;
    end
end

matrix = zeros(M*N,M*N);

%%rearrange the index of cell's coefficients
for i = 1 : M*N
    matrix(i,i)   =  re_aP(i) ;
end
for i = M+1 : M*N
    matrix(i,i-M) = -re_aW(i) ;
end
for i = 2 : M*N
    matrix(i,i-1) = -re_aN(i) ;
end
for i = 1 : M*N-1
    matrix(i,i+1) = -re_aS(i) ;
end
for i = 1 : M*N-M
    matrix(i,i+M) = -re_aE(i) ;
end

%draw
phi = matrix\re_Su;
%将phi再调回二维
for j = 1 : N
    for i = 1 : M
        re_phi(i,j) = phi((j-1)*M+i,1);
    end
end
figure(2);
clf;
mesh(x,y,re_phi');
colorbar;
xlabel('x(m)');
ylabel('y(m)');
zlabel('phi');
box on;
title(['2D FVM N-bc']);