function [u,x,y,t]=heat_ADI(a,D,T,ixy0,bxyt,Mx,My,N)
%a是系数，D是一组边界，T是总时间
%ixy0是初始条件(x,y,0)的映射关系，bxyt是(x，y，t)点的映射关系
%Mx，My，N分别是x方向，y方向，t方向上分成的份数

%首先 确定步长
%x
dx=D(1)/Mx;x=[0:Mx]*dx;
dy=D(2)/My;y=[0:My]'*dy;
dt=T/N;t=[0:N]*dt;
Mx1=Mx+1;
My1=My+1;

%initial condition
u=zeros(Mx-1,My-1);
for j=1:Mx-1
    for i=1:My-1
        u(i,j)=ixy0(x(i),y(j));
    end
end

%系数
rx=a*dt/dx/dx;
ry=a*dt/dy/dy;
Ax=zeros(Mx-1,My-1);
Ay=zeros(Mx-1,My-1);

%for x direction
for i=1:Mx-1
    Ax(i,i)=1+2*ry;
    if i>1
        Ax(i-1,i)=-ry;
        Ax(i,i-1)=-ry;
    end
end

%for y direction
for i=1:My-1
    Ay(i,i)=1+2*rx;
    if i>1
        Ay(i-1,i)=-rx;
        Ay(i,i-1)=-rx;
    end
end

%boundary condition
for k=1:N
    %u0=u;
    for j=1:My1
        u([1 Mx1],j)=[feval(bxyt,x(1),y(j),t(k)) feval(bxyt,x(Mx1),y(j),t(k))];
    end
    for i=1:Mx1
        u(i,[1 My1])=[feval(bxyt,x(i),y(1),t(k)) feval(bxyt,x(i),y(My1),t(k))];
    end
    u0=u;
    if mod(k,2)==0
        for i=2:My
            j=2:Mx;
            bx=[ry*u(i,1) zeros(1,Mx-3) ry*u(i,My1)]+rx*(u0(i-1,j)+u0(i+1,i))+(1-2*rx)*u0(i,j);
            u(i,j)=Ay\bx';
        end
    else
        for j=2:Mx
            i=2:My;
            by=[rx*u(1,j);zeros(My-3,1);rx*u(Mx1,j)]+ry*(u0(i,j-1)+u0(i,j+1))+(1-2*ry)*u0(i,j);
            u(i,j)=Ax\by;
        end
    end
    mesh(x,y,u);
    title(['step = ', num2str(k, '%d'), ', time = ', num2str(t(k), '%f')]);
    drawnow;
end