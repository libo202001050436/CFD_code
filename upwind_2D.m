function [u,aP,aE,aW,aN,aS,Su]=upwind_2D(rho,M,N,k,phix0,phixf,phiy0,phiyf,Lx,Ly)


dx = Lx  / M ;
dy = Ly  / N ;
Dx = k   / dx;
Dy = k   / dy;

x = zeros(M,1  );
y = zeros(1  ,N);

x(1:M-1,1) = dx/2 + [0:M-2] * dx ;
x(M  ,1) = Lx ;

y(1,1:N-1) = dy/2 + [0:N-2] * dy ;
y(1,  N) = Ly ;
u = zeros(M, N);
for i = 1:M
    u(i, [ 1 N ] ) = [phiy0(x(i)) phiyf(x(i))];
end
for j = 1:N
    u([1 M], j) = [phix0(y(j)) phixf(y(j))];
end
u(1, 1)     = (u(2, 1) + u(1, 2)) / 2;
u(1, N)     = (u(1, N-1) + u(2, N)) / 2;
u(M, 1)     = (u(M-1, 1) + u(M, 2)) / 2;
u(M, N)     = (u(M-1, N) + u(M, N-1)) / 2;
v = zeros(M,N);
v(1, 1) = (v(2, 1) + v(1, 2)) / 2;
v(1, N) = (v(1, N-1) + v(2, N)) /2;
v(M, 1) = (v(M-1, 1) + v(M, 2)) / 2;
v(M, N) = (v(M-1, N) + v(M, N-1)) /2;
Fx = rho * u ;
Fy = rho * v ;

%%
%aW = Dx + Fx ;
x = zeros(M, 1 );
y = zeros(1 , N);

x(1:M,1) = dx/2 + [0:M-1] * dx ;

y(1,1:N) = dy/2 + [0:N-1] * dy ;

aW = zeros(M,N) ;
% aW(1:M,1) = 0 ;
% aW(1:M,2:N) = Dx + Fx(1:M,2:N);
for i = 1:M
    for j = 1:N
        if j == 1
            aW(i,j) = 0;
        else
            aW(i,j) = Dx + max(Fx(i,j),0);
        end
    end
end

aE = zeros(M,N) ;
% aE(1:M,1:N-1) = Dx;
% aE(2:M, N)  = 0 ;
for i = 1:M
    for j = 1:N
        if j == N
            aE(i,j) = 0;
        else
            aE(i,j) = Dx + max(-Fx(i,j),0);
        end
    end
end

aN = zeros(M,N) ;
% aN(1,1:N) = 0;
% aN(2:M,1:N) = Dy -Fy(2:M,1:N);
for i = 1:M
    for j = 1:N
        if i == 1
            aN(i,j) = 0;
        else
            aN(i,j) = Dy + max(-Fy(i,j),0);
        end
    end
end

aS = zeros(M,N) ;
% aS(1:M-1,1:N) = Dy;
% aS(  M  ,1:N) = 0 ;
for i = 1:M
    for j = 1:N
        if i == M
            aS(i,j) = 0;
        else
            aS(i,j) = Dy + max(Fy(i,j),0);
        end
    end
end

SP = zeros(M,N) ;
SP(1:M,1) = -(2*Dx+Fx(1,1:M)) ;
SP(1:M,N) = -2*Dx ;
SP(1,1:N) = -(2*Dy+Fy(1,1:N)) ;
SP(M,1:N) = -2*Dy ;
SP( 1 , 1 ) = SP(1,2)+SP(2,1) ;
SP( 1 , N ) = SP(1,N-1)+SP(2,N);
SP( M , 1 ) = SP(M,2)+SP(M-1,1);
SP( M , N ) = SP(M-1,N)+SP(M,N-1);

aP = aW + aE + aS + aN -SP ;
clf;
mesh(x,y,aP);

Su = zeros(M,N) ;
for i = 2:N-1
    Su(1,i) = (2*Dy+Fy(i,j))*phiy0(x(1),y(i));
end
for i = 2:N-1
    Su(M,i) = (2*Dy)*phiyf(x(M),y(i));
end
for i = 2:M-1
    Su(i,1) = (2*Dx+Fx(i,j))*phix0(x(i),y(1));
end
for i = 2:M-1
    Su(i,N) = (2*Dx)*phixf(x(i),y(N));
end
Su( 1 , 1 ) = (Su(1,2)+Su(2,1));
Su( 1 , N ) = (Su(1,N-1)+Su(2,N));
Su( M , 1 ) = (Su(M-1,1)+Su(M,2));
Su( M , N ) = (Su(M-1,N)+Su(M,N-1));

%%
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

% matrix = zeros(M*N,M*N);

%%rearrange the index of cell's coefficients
% for i = 1 : M*N
%     matrix(i,i)   =  re_aP(i) ;
% end
% for i = M+1 : M*N
%     matrix(i,i-M) = -re_aW(i) ;
% end
% for i = 2 : M*N
%     matrix(i,i-1) = -re_aN(i) ;
% end
% for i = 1 : M*N-1
%     matrix(i,i+1) = -re_aS(i) ;
% end
% for i = 1 : M*N-M
%     matrix(i,i+M) = -re_aE(i) ;
% end
%sss = M*N + (M*N - M) + (M*N - 1) + (M * N - M) + (M*N-1)
p = 0;
for i = 1:M*N
    for n = 1:5
        if n ==1 && i > M+1
            p = p + 1;
            e(p) = -re_aW(i);
            r(p) = i ;
            c(p) = i-N ;
        end
        if n == 2 && i >= 2
            p = p + 1;
            e(p) = -re_aN(i);
            r(p) = i ;
            c(p) = i-1 ;
        end
        if n == 3
            p = p + 1;
            e(p) = re_aP(i);
            r(p) = i ;
            c(p) = i ;
        end
        if n == 4 && i <= M*N-1
            p = p + 1;
            e(p) = -re_aS(i);
            r(p) = i ;
            c(p) = 1+i ;
        end
        if n == 5 && i<M*N-M
            p = p + 1;
            e(p) = -re_aE(i);
            r(p) = i ;
            c(p) = M+i ;
        end
    end
end
matrix = sparse(r,c,e);
%full(matrix)

%draw
phi = matrix\re_Su;
%将phi再调回二维
re_phi = zeros(M,N);
% re_phi(2:M+1,1)=phix0;
% re_phi(2:M+1,N+2)=phixf;
% re_phi(1,2:N+1)=phiy0;
% re_phi(M+2,2:N+1)=phiyf;
% re_phi( 1 , 1 ) = (re_phi(1,2)+re_phi(2,1))/2;
% re_phi( 1 , N+2 ) = (re_phi(1,N+1)+re_phi(2,N+2))/2;
% re_phi( M+2 , 1 ) = (re_phi(M+1,1)+re_phi(M+2,2))/2;
% re_phi( M+2 , N+2 ) = (re_phi(M+1,N+2)+re_phi(M+2,N+1))/2;
for j = 1 : N
    for i = 1 : M
        re_phi(i,j) = phi((j-1)*M+i,1);
    end
end
figure(1);
clf;
mesh(x,y,re_phi');
xlabel('x(m)');
ylabel('y(m)');
zlabel('phi');
box on;
title('upwind 2D Dirichlet');
drawnow;
b=full(matrix);
figure(2);
c = M*N;
dx = 1/c;
x = [1:c]*dx;
mesh(x,x,b);
box on;
title('matrix');
drawnow;
u=phi;