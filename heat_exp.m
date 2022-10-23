function [u,x,t]=heat_exp(a,xf,T,it0,bx0,bxf,M,N)
%先将x和t分别插入M和N个间隔
dx=xf/M;x=[0:M]'*dx;
dt=T/N;t=[0:N]*dt;
%输入边界条件
M1=M+1;N1=N+1;
for m=1:M,u([1 M1],m)=[bx0(t(m)) bxf(t(m))];end
for n=1:N,u(n,1)=[it0(x(n))];end
dx2=dx*dx;r=a*dt/dx2;
%热传导方程用差分法差分后的公式
for tol=1:MaxIter
    for k=1:N
        for i=2:M
            u(i,k+1)=r*(u(i+1,k)+u(i-1,k))+(1-2*r)*u(i,k);
        end
    end
    if tol>max(max(abs(u-u0))),break;end
    u=u0;
end