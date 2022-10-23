function [u,x,t]=heat_imp(a,xf,T,it0,bx0,bxf,M,N)
%先将x和t分别插入M和N个间隔
dx=xf/M;x=[0:M]'*dx;
dt=T/N;t=[0:N]*dt;
%输入边界条件
M1=M+1;N1=N+1;
for m=1:M1,u([1 M1],m)=[bx0(t(m)) bxf(t(m))];end
for n=1:N1,u(n,1)=[it0(x(n))];end
dx2=dx*dx;
r=a*dt/dx2;
r2=1+2*r;
%热传导方程用差分法差分后的公式
for i=1:M-1
    A(i,i)=r2;
    if i>1,A(i-1,i)=-r;
        A(i,i-1)=-r;
    end
end

for k = 2:N1
    b = [r*u(1,k);zeros(M-3,1);r*u(M1,k)]+u(2:M,k-1);%材料上的下标很烂
    u(2:M,k)=A/b;%A是对角矩阵
end