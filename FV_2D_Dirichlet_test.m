clear;
clc;
Lx = 1 ;
Ly = 1 ;
Tx0 = 3 ;
Txf = 3 ;
Ty0 = 3 ;
Tyf = 3 ;
M = 100 ;
N = 110 ;
k = 1 ;
A = 1 ;
q = 1 ;
[u,x,y]=FV_2D_Dirichlet(Lx,Ly,Tx0,Txf,Ty0,Tyf,M,N,k,A,q);