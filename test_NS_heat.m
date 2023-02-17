clear;
clc;
Lx = 1;
Ly = 1;
time = 1;
dt = 0.01;
T0 = 0;
alpha = 0.01;
rho0 = 1;
dx = Lx/100;
dy = Ly/100;
mu = 1; %动力粘滞系数
nu = mu / rho0;%运动粘滞系数
T_char = 1;
L_char = 1;
V_char = L_char / T_char;
temperature_char = 1;
g = 1e7;
P_char = 1;
k = 1;
uN = @(x)(0);
uS = @(x)(0);
vW = @(y)(0);
vE = @(y)(0);
TS = @(x)(1);
TN = @(x)(0);
TW = @(y)(0);
TE = @(y)(0);
%WP means weighting parameter, if WP equal 1 means fully implicit, if WP
%equal 0 means explict, if WP equal 1/2 means Crank-Nicolson scheme
WP = 1/2;
heat_cap = 1;
%一定检查Re是否处于一个比较合理的范围
Re = L_char * L_char / T_char * nu;
%Ra大于10000是才有比较好的对流的效果
Ra = g * temperature_char * T_char * L_char * alpha / nu;
[T,u,v] = unsteady_heat(Lx,Ly,time,dt,T0,TW,TE,TN,TS,dx,dy,rho0,...
    alpha,V_char,L_char,T_char,mu,g,P_char,temperature_char,heat_cap,k,uN,uS,vW,vE,WP);