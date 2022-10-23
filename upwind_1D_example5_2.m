clear;
clc;

u   = 0.1 ;
Rho = 1   ;
k   = 0.1 ;
M   = 5   ;
L   = 1   ;
phi0= 1   ;
phiL= 0   ;
[U,x] = upwind_1D(Rho,M,k,u,phi0,phiL,L);
for i=1:1000+1
    dx=L/1000;
    x=[0:1000]*dx;
    phi(i) = (exp(Rho*u*x(i)/k)-1)/(exp(Rho*u*L/k)-1)*(phiL-phi0)+phi0;
end
plot(x,phi);
hold on;
A=1;
T0=phi0;
T00=phiL;
[u,x]=FV_1D_Neumam(L,T0,T00,M,n2);