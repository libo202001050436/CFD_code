function [u,x,y]=SIMPLE(uW,uE,uN,uS,vW,vE,vN,vS,...
    rho,pW,pE,pN,pS,alpha_p,alpha_u,Lx,Ly,M,N,mu)

dx = Lx / M;
dy = Ly / N;

x = zeros(M,1);
y = zeros(1,N);
x(1:M,1) = dx/2 + [0:M-1] * dx;

y(1,1:N) = dy/2 + [0:N-1] * dy;

%%boundary conduction
u = zeros(M, N);
% for i = 1:M
%     u(i, [ 1 N ] ) = [uW(x(i)) uE(x(i))];
% end
% for j = 1:N
%     u([1 M], j) = [uN(y(j)) uS(y(j))];
% end
u(1, 1)     = 0;
u(1, N)     = 0;
u(M, 1)     = 0;
u(M, N)     = 0;

v = zeros(M,N);
% for i = 1:M
%     v(i, [1 N]) = [vW(x(i)) vE(x(i))];
% end
% for j = 1:N
%     v([1 M], j) = [vN(y(j)) vS(y(j))];
% end
v(1, 1) = 0;
v(1, N) = 0;
v(M, 1) = 0;
v(M, N) = 0;

p = zeros(M,N);
for i = 1:M
    p(i, [1 N]) = [pW(x(i)) pE(x(i))];
end
for j = 1:N
    p([1 M], j) = [pN(y(j)) pS(y(j))];
end
p(1, 1) = (p(2, 1 ) + p(1, 2)) /2;
p(1, N) = (p(1, N-1) + p(2, N)) /2;
p(M, 1) = (p(M-1, 1) + p(M, 2)) /2;
p(M, N) = (p(M-1, N) + p(M, N-1)) /2;

u_cor = zeros(M,N);
v_cor = zeros(M,N);

ue = u;
uw = u;
vn = v;
vs = v;
b  = zeros(M,N);

norm_u0 = zeros(10000,1);
norm_v0 = zeros(10000,1);
norm_b  = zeros(10000,1);
% for i = 1:M
%     v(i, [1 N]) = [vW(x(i)) vE(x(i))];
% end
% for j = 1:N
%     v([1 M], j) = [vN(y(j)) vS(y(j))];
% end
% for i = 1:M
%     u(i, [1 N]) = [uW(x(i)) uE(x(i))];
% end
% for j = 1:N
%     u([1 M], j) = [uN(y(j)) uS(y(j))];
% end
%%开始迭代

for n = 0:3000

    for i = 1:M
        v(i, [1 N]) = [vW(x(i)) vE(x(i))];
    end
    for j = 1:N
        v([1 M], j) = [vN(y(j)) vS(y(j))];
    end
    for i = 1:M
        u(i, [1 N]) = [uW(x(i)) uE(x(i))];
    end
    for j = 1:N
        u([1 M], j) = [uN(y(j)) uS(y(j))];
    end
        u(1,1) = 0;
        u(1,N) = 0;
        v(1,1) = 0;
        v(1,N) = 0;

    p0 = p;
    u0 = u;
    v0 = v;
    p_cor = zeros(M*N,1);
    b     = zeros(M,N);

    pE_ghost(1:M) = p0(1:M, N) * 2 - p0(1:M, N-1);
    pW_ghost(1:M) = p0(1:M, 1) * 2 - p0(1:M, 2);
    pN_ghost(1:N) = p0(1, 1:N) * 2 - p0(2, 1:N);
    pS_ghost(1:N) = p0(M, 1:N) * 2 - p0(M-1, 1:N);
    uE_ghost(1:M,1) = u0(1:M, N) * 2 - u0(1:M, N-1);
    uW_ghost(1:M,1) = u0(1:M, 1) * 2 - u0(1:M, 2);
    uN_ghost(1,1:N) = u0(1, 1:N) * 2 - u0(2, 1:N);
    uS_ghost(1,1:N) = u0(M, 1:N) * 2 - u0(M-1, 1:N);
    vE_ghost(1:M,1) = v0(1:M, N) * 2 - v0(1:M, N-1);
    vW_ghost(1:M,1) = v0(1:M, 1) * 2 - v0(1:M, 2);
    vN_ghost(1,1:N) = v0(1, 1:N) * 2 - v0(2, 1:N);
    vS_ghost(1,1:N) = v0(M, 1:N) * 2 - v0(M-1, 1:N);
    %     uE_ghost(1:M,1) = uE(x(1:M));
    %     uW_ghost(1:M,1) = uW(x(1:M));
    %     uN_ghost(1,1:N) = uN(y(1:N));
    %     uS_ghost(1,1:N) = uS(y(1:N));
    %     vE_ghost(1:M,1) = vE(x(1:M));
    %     vW_ghost(1:M,1) = vW(x(1:M));
    %     vN_ghost(1,1:N) = vN(y(1:N));
    %     vS_ghost(1,1:N) = vS(y(1:N));

    %%u0,v0为上一层迭代的数据,根据u0,v0确定动量u_star,v_star
    %先算u_star
    uu = u0;
    vv = v0;

    [u0,aPu,aEu,aWu,aNu,aSu,Suu]  = ...
        update_u_momentum(Lx,Ly,uu,ue,uw,vv,vn,vs,M,N,mu,rho,p0,alpha_u);

    %计算u_wave(P,E,W,N,S)
    u_wave = zeros(M,N);
    u_wave0= zeros(M,N);

    for i =1 : M
        for j = 1:N
            if i>1 && i<M && j>1 && j<N
                u_wave(i,j) = aEu(i,j) * u0(i,j+1) + aWu(i,j) *u0(i,j-1) + aSu(i,j) *...
                    u0(i+1,j) + aNu(i,j) * u0(i-1,j)+Suu(i,j);
                u_wave0(i,j) = u_wave(i,j) / aPu(i,j);
            end
            if i == 1 && j ~= 1 && j ~= N
                u_wave(i,j) = aEu(i,j) * u0(i,j+1) + aWu(i,j) *u0(i,j-1) + aSu(i,j) *...
                    u0(i+1,j) + uN_ghost(j) * aNu(i,j)+Suu(i,j);
                u_wave0(i,j) = u_wave(i,j) / aPu(i,j);
            end
            if i == M && j ~= 1 && j ~= N
                u_wave(i,j) = aEu(i,j) * u0(i,j+1) + aWu(i,j) *u0(i,j-1) + aSu(i,j) *...
                    uS_ghost(j) + aNu(i,j) * u0(i-1,j)+Suu(i,j);
                u_wave0(i,j) = u_wave(i,j) / aPu(i,j);
            end
            if j == 1 && i ~= 1 && i ~= M
                u_wave(i,j) = aEu(i,j) * u0(i,j+1) + uW_ghost(i) *aWu(i,j) + aSu(i,j) *...
                    u0(i+1,j) + aNu(i,j) * u0(i-1,j)+Suu(i,j);
                u_wave0(i,j) = u_wave(i,j) / aPu(i,j);
            end
            if j == N && i ~= 1 && i ~= M
                u_wave(i,j) = uE_ghost(i) * aEu(i,j) + aWu(i,j) *u0(i,j-1) + aSu(i,j) *...
                    u0(i+1,j) + aNu(i,j) * u0(i-1,j)+Suu(i,j);
                u_wave0(i,j) = u_wave(i,j) / aPu(i,j);
            end
            if j == 1 && i == 1
                u_wave(1,1) = uW_ghost(1)*aWu(1,1) + aEu(1,1)*u0(i,j+1) + aNu(i,j)*uN_ghost(j)...
                    +aSu(i,j)*uS_ghost(1)+Suu(i,j);
                u_wave0(i,j) = u_wave(i,j) / aPu(i,j);
            end
            if i == M && j==1
                u_wave(i,j) = aWu(i,j)*uW_ghost(i) + aEu(i,j)*u0(i,j+1) + aNu(i,j)*u0(i-1,j)...
                    +aSu(i,j)*uS_ghost(j)+Suu(i,j);
                u_wave0(i,j) = u_wave(i,j) / aPu(i,j);
            end
            if i == 1 && j ==N
                u_wave(i,j) = aWu(i,j)*u0(i,j-1) + aEu(i,j)*uE_ghost(i) + aNu(i,j)*uN_ghost(j)...
                    +aSu(i,j)*uS_ghost(j)+Suu(i,j);
                u_wave0(i,j) = u_wave(i,j) / aPu(i,j);
            end
            if i == M && j ==N
                u_wave(i,j) = aWu(i,j)*u0(i,j-1) + aEu(i,j)*uE_ghost(i) + aNu(i,j)*u0(i-1,j)...
                    +aSu(i,j)*uS_ghost(j)+Suu(i,j);
                u_wave0(i,j) = u_wave(i,j) / aPu(i,j);
            end
        end
    end
    u_wave = u_wave0;

    %     for i = 1:M
    %         for j = 1:N
    %             if j > 1 && j < N
    %                 u0(i,j) = u_wave(i,j) - dy / aPu(i,j) * (p0(i,j+1) - p0(i,j-1)) / 2;
    %             end
    %             if j == 1
    %                 u0(i,1) = u_wave(i,1) - dy / aPu(i,1) * (p0(i,j+1) - pW_ghost(j)) / 2;
    %             end
    %             if j == N
    %                 u0(i,N) = u_wave(i,j) - dy / aPu(i,j) * (pE_ghost(j) - p0(i,j-1)) / 2;
    %             end
    %         end
    %     end
    uW_wave_ghost(1:M,1) = u_wave(1:M, 1) * 2 - u_wave(1:M, 2);
    uE_wave_ghost(1:M,1) = u_wave(1:M, N) * 2 - u_wave(1:M, N-1);

    %计算v_star

    [v0,aPv,aEv,aWv,aNv,aSv,Suv]  = ...
        update_v_momentum(Ly,Lx,uu,ue,uw,vv,vn,vs,M,N,mu,rho,p0,alpha_u);

    %计算v_wave(P,E,W,N,S)
    v_wave = zeros(M,N);

    %     v_wave(1,1) = aEv(1,1) * v0(1,2) + aWv(1,1) *vW_ghost(1) + aSv(1,1) *...
    %         v0(2,1) + aNv(1,1) * vN_ghost(1);
    %     v_wave(1,N) = aEv(1,N) * vE_ghost(N) + aWv(1,N) *v0(1,N) + aSv(1,N) *...
    %         v0(2,N) + aNv(1,N) * vN_ghost(N);
    %     v_wave(M,1) = aEv(M,1) * vE_ghost(M) + aWv(M,1) *v0(1,N) + aSv(M,1) *...
    %         vS_ghost(1) + aNv(M,1) * v0(M-1,1);
    %     v_wave(M,N) = aEv(M,N) * vE_ghost(N) + aWv(M,N) *v0(M,N) + aSv(M,N) *...
    %         vS_ghost(N) + aNv(M,N) * v0(M-1,N);
    v_wave0=zeros(M,N);
    for i =1 : M
        for j = 1:N
            if i>1 && i<M && j>1 && j<N
                v_wave(i,j) = aEv(i,j) * v0(i,j+1) + aWv(i,j) *v0(i,j-1) + aSv(i,j) *...
                    v0(i+1,j) + aNv(i,j) * v0(i-1,j)+Suv(i,j);
                v_wave0(i,j) = v_wave(i,j)/aPv(i,j);
            end
            if i == 1 && j > 1 && j < N
                v_wave(i,j) = aEv(i,j) * v0(i,j+1) + aWv(i,j) *v0(i,j-1) + aSv(i,j) *...
                    v0(i+1,j) + vN_ghost(j) * aNv(i,j)+Suv(i,j);
                v_wave0(i,j) = v_wave(i,j)/aPv(i,j);
            end
            if i == M && j > 1 && j < N
                v_wave(i,j) = aEv(i,j) * v0(i,j+1) + aWv(i,j) *v0(i,j-1) + aSv(i,j) *...
                    vS_ghost(j) + aNv(i,j) * v0(i-1,j)+Suv(i,j);
                v_wave0(i,j) = v_wave(i,j)/aPv(i,j);
            end
            if j == 1 && i > 1 && i < M
                v_wave(i,j) = aEv(i,j) * v0(i,j+1) + vW_ghost(i) *aWv(i,j) + aSv(i,j) *...
                    v0(i+1,j) + aNv(i,j) * v0(i-1,j)+Suv(i,j);
                v_wave0(i,j) = v_wave(i,j)/aPv(i,j);
            end
            if j == N && i > 1 && i < M
                v_wave(i,j) = vE_ghost(i) * aEv(i,j) + aWv(i,j) *v0(i,j-1) + aSv(i,j) *...
                    v0(i+1,j) + aNv(i,j) * v0(i-1,j)+Suv(i,j);
                v_wave0(i,j) = v_wave(i,j)/aPv(i,j);
            end
            if i == 1 && j == 1
                v_wave(i,j) = aEv(i,j)*v0(i,j+1) + aWv(i,j)*vW_ghost(i) + aNv(i,j)*vN_ghost(j)...
                    +aSv(i,j)*v0(i+1,j)+Suv(i,j);
                v_wave0(i,j) = v_wave(i,j)/aPv(i,j);
            end
            if i == 1 && j == N
                v_wave(i,j) = aEv(i,j)*vE_ghost(i) + aWv(i,j)*v0(i,j-1) + aNv(i,j)*vN_ghost(j)...
                    +aSv(i,j)*v0(i+1,j)+Suv(i,j);
                v_wave0(i,j) = v_wave(i,j)/aPv(i,j);
            end
            if i == M && j == 1
                v_wave(i,j) = aEv(i,j)*v0(i,j+1) + aWv(i,j)*vW_ghost(i) + aNv(i,j)*v0(i-1,j)...
                    +aSv(i,j)*vS_ghost(j)+Suv(i,j);
                v_wave0(i,j) = v_wave(i,j)/aPv(i,j);
            end
            if i == M && j == N
                v_wave(i,j) = aEv(i,j)*vE_ghost(i) + aWv(i,j)*v0(i,j-1) + aNv(i,j)*v0(i-1,j)...
                    +aSv(i,j)*vS_ghost(j)+Suv(i,j);
                v_wave0(i,j) = v_wave(i,j)/aPv(i,j);
            end
        end
    end
    v_wave=v_wave0;
    %     for i = 1:M
    %         for j = 1:N
    %             if i > 1 && i < M
    %                 v0(i,j) = v_wave(i,j) - dx / aPv(i,j) * (p0(i-1,j) - p0(i+1,j)) / 2;
    %             end
    %             if i == 1
    %                 v0(i,1) = v_wave(i,j) - dx / aPv(i,j) * (pN_ghost(j) - p0(i+1,j)) / 2;
    %             end
    %             if i == N
    %                 v0(i,N) = v_wave(i,j) - dx / aPv(i,j) * (p0(i-1,j) - pS_ghost(j))/2;
    %             end
    %         end
    %     end
    vN_wave_ghost(1:N) = v_wave(1, 1:N) * 2 - v_wave(2, 1:N);
    vS_wave_ghost(1:N) = v_wave(M, 1:N) * 2 - v_wave(M-1, 1:N);

    %%计算(ue，uw，vn，vs)


    for i = 1:M
        for j = 1:N
            if j < N
                ue(i,j) = (u_wave(i,j)+u_wave(i,j+1))/2 - 2*dy / (aPu(i,j)+aPu(i,j+1)) * (p0(i,j+1) - p0(i,j));
            end
            if j == N
                aPu_ghost = 2*aPu(i,j) - aPu(i,j-1);
                ue(i,j) = (u_wave(i,j)+uE_wave_ghost(i))/2 - 2*dy / (aPu(i,j)+aPu_ghost) * (pE_ghost(i) - p0(i,j));
            end
            if j > 1
                uw(i,j) = (u_wave(i,j)+u_wave(i,j-1))/2 - 2*dy / (aPu(i,j)+aPu(i,j-1)) * (p0(i,j) - p0(i,j-1));
            end
            if j ==1
                aPu_ghost = 2*aPu(i,j) - aPu(i,j+1);
                uw(i,j) = (u_wave(i,j)+uW_wave_ghost(i))/2 - 2*dy / (aPu(i,j)+aPu_ghost) * (p0(i,j) - pW_ghost(i));
            end
            if i > 1
                vn(i,j) = (v_wave(i,j)+v_wave(i-1,j))/2 - 2*dx / (aPv(i,j)+aPv(i-1,j)) * (p0(i-1,j) - p0(i,j));
            end
            if i == 1
                aPv_ghost = 2*aPv(i,j) - aPv(i+1,j);
                vn(i,j) = (v_wave(i,j)+vN_wave_ghost(j))/2 - 2*dx / (aPv(i,j)+aPv_ghost) * (pN_ghost(j) - p0(i,j));
            end
            if i < M
                vs(i,j) = (v_wave(i,j)+v_wave(i+1,j))/2 - 2*dx / (aPv(i,j)+aPv(i+1,j)) * (p0(i,j) - p0(i+1,j));
            end
            if i == M
                aPv_ghost = 2*aPv(i,j) - aPv(i-1,j);
                vs(i,j) = (v_wave(i,j)+vS_wave_ghost(1,j))/2 - 2*dx / (aPv(i,j)+aPv_ghost) * (p0(i,j) - pS_ghost(j));
            end
        end
    end
    %%
    ue(1:M,N) = zeros(M,1);
    uw(1:M,1) = zeros(M,1);
    vn(1,1:N) = zeros(1,N);
    vs(M,1:N) = zeros(1,N);

    %压力修正值源项
    %         for i = 1:M
    %             for j = 1:N
    %                 if i >1 && i <M && j >1 && j < N
    %                     b(i,j) = rho * (uw(i,j-1) * dy - ue(i,j+1) * dy + vs(i-1,j) * dx - vn(i+1,j) * dx);
    %                 end
    %                 if i == 1 && j>1 && j<N
    %                     b_ghost = vn(1,j) * 2 - vn(2,j);
    %                     b(i,j) = rho * (vs(i+1,j)  * dx - b_ghost * dx + uw(i,j-1) * dy - ue(i,j+1) * dy);
    %                 end
    %                 if i == M && j>1 && j<N
    %                     b_ghost = vs(M,j) * 2 - vs(M-1,j);
    %                     b(i,j) = rho * (uw(i,j-1) * dy - ue(i,j+1) * dy - vn(i-1,j) * dx + b_ghost * dx);
    %                 end
    %                 if j == 1 && i >1 && i <M
    %                     b_ghost = uw(i,j)*2 - uw(i,j+1);
    %                     b(i,j) = rho * (b_ghost * dy - ue(i,j+1) * dy + vs(i-1,j) * dx - vn(i+1,j) * dx );
    %                 end
    %                 if j == N && i >1 && i <M
    %                     b_ghost = ue(i,N) * 2 - ue(i,N-1);
    %                     b(i,j) = rho * (uw(i,j-1) * dy - b_ghost * dy + vs(i-1,j) * dx - vn(i+1,j) * dx);
    %                 end
    %                 if i == 1 && j == 1
    %                     b_ghost = -vn(i,j) * 2 * dx - uw(i,j+1) * dy + vn(i+1,j) * dx + uw(i,j) * 2 * dy;
    %                     b(i,j) = rho * (b_ghost - ue(i,j+1) * dy - vn(i+1,j) * dx);
    %                 end
    %                 if i == M && j == N
    %                     b_ghost = -vs(i,j) * 2 * dx + ue(i,j) * dy * 2 +vs(i-1,j) * dx - ue(i,j-1) * dy;
    %                     b(i,j) = rho * (uw(i,j-1) * dy - vn(i-1,j) * dx - b_ghost);
    %                 end
    %                 if i == 1 && j == N
    %                     b_ghost = ue(i,j) * 2 * dy - ue(i,j-1) * dy + vn(i,j) * 2 * dx - vn(i+1,j) * dx;
    %                     b(i,j) = rho * (uw(i,j-1) * dy - b_ghost + vs(i+1,j) * dx);
    %                 end
    %                 if i == M && j ==1
    %                     b_ghost = uw(i,j) * 2 * dy - uw(i,j+1) * dy + vs(i,j) * 2 * dx - vs(i,j+1) * dx;
    %                     b(i,j) = rho * (b_ghost - ue(i,j+1) * dy - vn(i-1,j) * dx);
    %                 end
    %             end
    %         end
    for i = 1:M
        for j = 1:N
            b(i,j) = rho * (uw(i,j) - ue(i,j) + vs(i,j) - vn(i,j));
        end
    end
    %压力修正值p_cor
    de = zeros(M,N);
    ds = zeros(M,N);
    dn = zeros(M,N);
    dw = zeros(M,N);
    for i = 1:M
        for j = 1:N
            if i >1
                dn(i,j) = dy * 0.5 * (1/aPv(i,j) + 1 / aPv(i-1,j));
            end
            if i == 1
                dn(i,j) = dy * 0.5 * ( 1/aPv(i,j) + 1/(aPv(i,j) * 2 - aPv(i+1,j)));
            end
            if i < M
                ds(i,j) = dy * 0.5 * (1/aPv(i,j) + 1 / aPv(i+1,j));
            end
            if i == M
                ds(i,j) = dy * 0.5 * ( 1/aPv(i,j) + 1/(aPv(i,j) * 2 - aPv(i-1,j)));
            end
            if j < M
                de(i,j) = dx * 0.5 * (1/aPu(i,j) + 1 / aPu(i,j+1));
            end
            if j == M
                de(i,j) = dx * 0.5 * (1/aPu(i,j) + 1/(aPu(i,j) * 2 - aPu(i,j-1)));
            end
            if j > 1
                dw(i,j) = dx * 0.5 * (1/aPu(i,j) + 1 / aPu(i,j-1));
            end
            if j == 1
                dw(i,j) = dx * 0.5 * (1/aPu(i,j) + 1 / (aPu(i,j) * 2 - aPu(i,j+1)));
            end
            %             if i == 1
            %                 aP_ghost = aPv(i,j) * 2 -aPv(i+1,j);
            %                 dn(i,j) = 0.5 * dx * (1/aPv(i,j) + 1 / aP_ghost);
            %             end
            %             if i == M
            %                 aP_ghost = aPv(i,j) * 2 - aPv(i-1,j);
            %                 ds(i,j) = 0.5 * dx * (1/aPv(i,j) + 1 / aP_ghost);
            %             end
            %             if j == 1
            %                 aP_ghost = aPu(i,j) * 2 - aPu(i,j+1);
            %                 dw(i,j) = 0.5 * dy * (1/aPu(i,j) + 1 / aP_ghost);
            %             end
            %             if j == N
            %                 aP_ghost = aPu(i,j) * 2 - aPu(i,j-1);
            %                 de(i,j) = 0.5 * dy * (1/aPu(i,j) + 1 / aP_ghost);
            %             end
        end
    end

    aE = rho * de ;
    aW = rho * dw ;
    aN = rho * dn ;
    aS = rho * ds ;
    aP = aE + aN + aS + aW;
    aE(1:M,N) = zeros(M,1);
    aN(1,1:N) = zeros(1,N);
    aS(M,1:N) = zeros(1,N);
    aW(1:M,1) = zeros(M,1);
    %     aP(1,N) = 1e10;
%         aP(1,1) = 1e10;
%         aP(M,1) = 1e10;
%         aP(M,N) = 1e10;
    [matrix] = sparse_coef_auto(aP,aW,aE,aN,aS,M,N);
    back_p_cor = zeros(M,N);
    b = reshape_a(b,M,N);
    p_cor = matrix \ b;
    for j = 1 : N
        for i = 1 : M
            back_p_cor(i,j) = p_cor((j-1)*M+i,1);%dimension up
        end
    end
    back_b = zeros(M,N);
    for j = 1 : N
        for i = 1 : M
            back_b(i,j) = b((j-1)*M+i,1);%dimension up
        end
    end
    b=back_b;
    p_cor = back_p_cor;
    %计算速度修正值u_cor，和v_cor
    %     p = p0 + alpha * p_cor;


    for i = 1:M
        for j = 1:N
            if j < N && j > 1
                u_cor(i,j) = dx / aPu(i,j) * (p_cor(i,j-1) - p_cor(i,j+1))/2;
            end
            if j == N
                p_ghost = p_cor(i,j) * 2 - p_cor(i,j-1);
                %                 p_ghost = 0;
                u_cor(i,j) = dx / aPu(i,j) * (p_cor(i,j-1) - p_ghost)/2;
            end
            if j == 1
                p_ghost = p_cor(i,j) * 2 - p_cor(i,j+1);
                %                 p_ghost = 0;
                u_cor(i,j) = dx / aPu(i,j) * (p_ghost - p_cor(i,j+1))/2;
            end
            if i > 1 && i<M
                v_cor(i,j) = dy / aPv(i,j) * (p_cor(i+1,j) - p_cor(i-1,j))/2;
            end
            if i == 1
                p_ghost = p_cor(i,j) * 2 - p_cor(i+1,j);
                %                 p_ghost = 0;
                v_cor(i,j) = dy / aPv(i,j) * (p_cor(i+1,j) - p_ghost)/2;
            end
            if i == M
                p_ghost = p_cor(i,j) * 2 - p_cor(i-1,j);
                %                 p_ghost = 0;
                v_cor(i,j) = dy / aPv(i,j) * (p_ghost - p_cor(i-1,j))/2;
            end
        end
    end

    %     for i = 1:M
    %         for j = 1:N
    %             if j < N && j > 1
    %                 u_cor(i,j) = A / aPu(i,j) * (p(i,j-1) - p(i,j+1))/2;
    %             end
    %             if j == N
    %                 p_ghost = p(i,j) * 2 - p(i,j-1);
    %                 u_cor(i,j) = A / aPu(i,j) * (p(i,j-1) - p_ghost)/2;
    %             end
    %             if j == 1
    %                 p_ghost = p(i,j) * 2 - p(i,j+1);
    %                 u_cor(i,j) = A / aPu(i,j) * (p_ghost - p(i,j+1))/2;
    %             end
    %             if i > 1 && i<M
    %                 v_cor(i,j) = A / aPv(i,j) * (p(i+1,j) - p(i-1,j))/2;
    %             end
    %             if i == 1
    %                 p_ghost = p(i,j) * 2 - p(i+1,j);
    %                 v_cor(i,j) = A / aPv(i,j) * (p(i+1,j) - p_ghost)/2;
    %             end
    %             if i == M
    %                 p_ghost = p(i,j) * 2 - p(i-1,j);
    %                 v_cor(i,j) = A / aPv(i,j) * (p_ghost - p(i-1,j))/2;
    %             end
    %         end
    %    end
    %     v_cor(1:M,[1 N]) = [zeros(M,1) zeros(M,1)];
    %     v_cor([1 M],1:N) = [zeros(1,N);zeros(1,N)];
    p = p0 + alpha_p * p_cor;
    u = u0 + u_cor;
    v = v0 + v_cor;
    n = n + 1;
    norm_b(n)  = norm(b,2);
    norm_v0(n) = norm(v_cor,2);
    norm_u0(n) = norm(u_cor,2);
    k=n;
    %     %%
    % plot(n,u_cor);
    if mod(k,100) == 0
        figure(1);
        subplot(4,4,1);
        mesh(y,x,u0);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('u0');
        box on;
        title('u0');
        subplot(4,4,2);
        mesh(y,x,v0);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('v0');
        box on;
        title('v0');
        subplot(4,4,3);
        mesh(y,x,p0);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('p');
        box on;
        title('p0');
        subplot(4,4,5);
        mesh(y,x,u_wave);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('u wave');
        box on;
        title('u wave');
        subplot(4,4,6);
        mesh(y,x,v_wave);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('v wave');
        box on;
        title('v wave');
        subplot(4,4,9);
        mesh(y,x,ue);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('ue');
        box on;
        title('ue');
        subplot(4,4,4);
        matrix = full(matrix);
        dxx = 1/M/N;
        xx  = [0:M*N-1] * dxx;
        mesh(xx,xx,matrix);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('mtrx');
        box on;
        title('matrix');
        subplot(4,4,10);
        mesh(x,y,uw);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('uw');
        box on;
        title('uw');
        subplot(4,4,11);
        mesh(x,y,vn);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('vn');
        title('vn');
        subplot(4,4,12);
        mesh(x,y,vs);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('vn');
        title('vs');
        subplot(4,4,13);
        mesh(x,y,u_cor);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('u cor');
        title('u cor');
        subplot(4,4,14);
        mesh(x,y,v_cor);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('v cor');
        title('v cor');
        subplot(4,4,15);
        mesh(x,y,p_cor);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('p cor');
        title('p cor');
        subplot(4,4,16);
        mesh(x,y,b);
        xlabel('x(m)');
        ylabel('y(m)');
        zlabel('b');
        title('b');
        subplot(4,4,7);
        mesh(y,x,aPu);
        xlabel('x');
        ylabel('y');
        title('aPu');
        subplot(4,4,8);
        mesh(y,x,aPv);
        xlabel('x');
        ylabel('y');
        title('aPv');
        figure(2);
        plot(1:n,norm_u0(1:n),'g',1:n,norm_v0(1:n),'r',1:n,norm_b(1:n),'b');
        %     semilogy(n,10.^(-n));
        legend('u0','v0','b');
    end


    for i = 1:M
        u(i,1) = uW(y(i));
        u(i,N) = uE(y(i));
    end
    for i = 1:N
        u(1,i) = uN(x(i));
        u(M,i) = uS(x(i));
        v(1,i) = vN(x(i));
        v(M,i) = vS(x(i));
    end
    if norm_b(n)<= 1e-10
        break;
    end



end
for i = 1:M
    for j = 1:N
        u0(i,j) = u(M-i+1,j);
        v0(i,j) = v(M-i+1,j);
        p0(i,j) = p(M-i+1,j);
    end
end
figure(3);
clf;
pcolor(x,y,p0);
title('p');
% hold on;
figure(4);
quiver(x,y,u0,v0);
startx = linspace(0,1,15);
starty = startx;
streamline(x,y,u0,v0,starty,starty);
end