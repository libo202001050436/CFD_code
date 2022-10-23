clear;
clc;
rho   = 1;
M     = 10;
N     = 10;
k     = 0.01;
u     = 1;
v     = 0;
phix0 = @(x,y)(1);
phixf = @(x,y)(0);
phiy0 = @(x,y)(0);
phiyf = @(x,y)(0);
Lx    = 1;
Ly    = 1;

[phi,x,y]=upwind_2D(rho,M,N,k,phix0,phixf,phiy0,phiyf,Lx,Ly);