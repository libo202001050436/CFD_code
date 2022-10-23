%solve_example_9_3
clear;
clc;
clf;
a=1e-4;
ixy0=inline('0','x','y');
bxyt=inline('exp(y)*cos(x)-exp(x)*cos(y)','x','y','t');
D=[4 4];
T=5000;
Mx=40;
My=40;
N=50;
[u,x,y,t]=heat_ADI(a,D,T,ixy0,bxyt,Mx,My,N);
mesh(x,y,u)