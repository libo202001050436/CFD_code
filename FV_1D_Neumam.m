%FVM_1D_steady_state_diffusion
function [u,x]=FV_1D_Neumam(L,T0,T00,M,n2)
%A:cross-sectional area
%k:thermal conductivity
%M:divide length of fin L into M equal control volumes
%N:divide Temperature T into N
%S:source term
%q:heat generation
%h:heat transfer coefficient
%T00:ambient temperature

dx = L/M;
x  = dx/2+[0:M-1]*dx;

%aW
aW        = zeros(M,1);
aW(2:M,1) = 1/dx;

%aE
aE        = zeros(M,1);
aE(1:M-1) = 1/dx;

%SP
SP        =  zeros(M,1);
SP(1,1)   = -n2*dx-2/dx;
SP(2:M,1) = -n2*dx;

%aP
aP      = zeros(M,1);
for i=1:M
    aP(i,1) = aW(i,1)+aE(i,1)-SP(i,1);
end

%Su
Su        = zeros(M,1);
Su(1,1)   = n2*dx*T00+2/dx*T0;
Su(2:M,1) = n2*dx*T00;

%coefficient
%aW
for i = 1:M-1
    u(i+1,i) = -aW(i+1);
end
%aP
for i = 1:M
    u(i,i)   =  aP(i);
end
%aE
for i = 1:M-1
    u(i,i+1) = -aE(i);
end
T=u\Su;
plot(x,T);