function [U,x]=upwind_1D(Rho,M,k,u,phi0,phiL,L)
%Rho:density

F=Rho*u;
dx=L/M;
D=k/dx;

x = zeros(1,M+2);
x(1,2:M+1)=dx/2+[0:M-1]*dx;
x(1,M+2)=L;

aW = zeros(M,1);
aW(1,1) = 0 ;
aW([2:M],1) = D + F ;

aE = zeros(M,1);
aE(1:M-1,1) = D ;
aE(M,1) = 0;

SP = zeros(M,1);
SP(1,1) = -(2 * D + F );
SP(2:M-1,1) = 0;
SP(M,1) = -2*D;

aP = zeros(M,1);
aP = aW + aE - SP;

Su = zeros(M,1);
Su(1,1) = (2*D+F ) * phi0;
Su(2:M-1,1) = 0;
Su(M,1) = 2*D*phiL;

U=zeros( M );
%aP
for i = 1:M
    U(i,i) = aP(i);
end
%aW
for i = 2:M
    U(i,i-1) = -aW(i);
end
%aE
for i = 1:M-1
    U(i,i+1) = -aE(i);
end
phi = zeros(M+2,1);
phi(1) = phi0;
phi(M+2,1) = phiL;
phi(2:M+1,1) = U\Su;
plot(x,phi);
hold on;