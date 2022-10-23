%FVM_1D_steady_state_diffusion
function [u,x]=FV_1D_Dirichlet(A,k,L,T0,Tf,M,q)
%A:cross-sectional area
%k:thermal conductivity
%M:divide length of rod L into M equal control volumes
%N:divide Temperature T into N
%S:source term
%q:heat generation

%coefficient
dx=L/M;
x=dx/2+[0:M-1]*dx;

%a_p*phi_p=a_w*phi_w+a_E*phi_E+S_u
%a_w
u(2:M,1)=-k*A/dx;

%a_E
u(1:M-1,2)=-k*A/dx;

%S=S_u+S_p*phi
%S_p
u(2:M-1,4)=0;
u([1 M],4)=-2*k/dx;

%a_p=a_w+a_E-S_p
u(1:M,3)=-u(1:M,1)-u([1:M],2)-u(1:M,4);

%S_u
u(1,5)   = q * A * dx+2*k*A*T0/dx;
u(2:M,5)= q * A * dx;
u(M,5)  = q * A * dx+2*k*A*Tf/dx;

c=zeros(M,M);
%a_W
for i=1:M-1
    c(i+1,i)=u(i+1,1);
end
%a_p
for i=1:M
    c(i,i)=u(i,3);
end
%a_E
for i=2:M
    c(i-1,i)=u(i-1,2);
end

T=c\u(1:M,5);
plot(x,T);
hold on;
