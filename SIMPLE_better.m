function [u,v,p]=SIMPLE_better(Lx,Ly,M,N,uN,uS,vW,vE,Fr,Eu,Re,rho)

alpha = 0.6;
alpha_p = 0.5;

dx = Lx / M;
dy = Ly / N;

x = zeros(M,1);
y = zeros(1,N);

x(1:M,1) = dx/2 + [0:M-1] * dx;
y(1,1:N) = dy/2 + [0:N-1] * dy;

u = zeros(M,N);
v = zeros(M,N);
p = zeros(M,N);
b = zeros(M,N);
u_wave = zeros(M,N);
v_wave = zeros(M,N);
u_cor = zeros(M,N);
v_cor = zeros(M,N);
norm_u0 = zeros(10000,1);
norm_v0 = zeros(10000,1);
norm_b  = zeros(10000,1);



%%iter
for n = 0:1e10

%     u(1:M,1) = zeros(M,1);
%     u(1:M,N) = zeros(M,1);
%     v(1,1:N) = zeros(1,N);
%     v(M,1:N) = zeros(1,N);
    %     u(7,1) = 0.1;
    if n == 0
        ue = u;
        uw = u;
        vn = v;
        vs = v;
    end
%         %         for i = 1:M
%         %             for j = 1:N
%         %                 rho(i,j) = rho0;
%         %             end
%         %         end
%         if u(1,1) ~= 0 || v(1,1) ~= 0
%             vn(1,1) = -u(1,1)+v(1,1);
%             vs(1,1) = 0;
%             uw(1,1) = 0;
%             ue(1,1) = u(1,1)-v(1,1);
%         end
%         if u(1,N) ~= 0 || v(1,N) ~= 0
%             vn(1,N) = u(1,N)+v(1,N);
%             vs(1,N) = 0;
%             uw(1,N) = u(1,N)+v(1,N);
%             ue(1,N) = 0;
%         end
%         if u(M,1) ~= 0 || v(M,1) ~= 0
%             vn(M,1) = 0;
%             vs(M,1) = u(M,1)+v(M,1);
%             uw(M,1) = 0;
%             ue(M,1) = u(M,1)-v(M,1);
%         end
%         if u(M,N) ~= 0 || v(M,N) ~= 0
%             vn(M,N) = 0;
%             vs(M,N) = -u(M,N)+v(M,N);
%             uw(M,N) = u(M,N)-v(M,N);
%             ue(M,N) = 0;
%         end


    %因为边界条件的存在，所以ue，uw，vn，vs均要符合边界条件

for j = 1:N
    for i = 1:M
        if i == M
            u(i,j) = uN(y(j));
        end
        if i == 1
            u(i,j) = uS(y(j));
        end
        if j == 1
            v(i,j) = vW(x(i));
        end
        if j == N
            v(i,j) = vE(x(i));
        end
    end
end


    p0 = p;
    u0 = u;
    v0 = v;

    %calculate new u0,v0
    uu = u0;
    vv = v0;
    [u0,aPu,aEu,aWu,aNu,aSu,Suu] = update_momentum(Lx,Ly,...
        uu,ue,uw,vv,vn,vs,M,N,p,alpha,Fr,Eu,Re,1,rho);
    %1 means function will return u
    [v0,aPv,aEv,aWv,aNv,aSv,Suv] = update_momentum(Lx,Ly,...
        uu,ue,uw,vv,vn,vs,M,N,p,alpha,Fr,Eu,Re,2,rho);
    %2 means function will return v

    %pseudo u,v.here is u_wave,v_wave
    for j = 1:N
        for i = 1:M
            if i == 1
                if j == 1
                    uu_ghost = u0(i,j) * 2 - u0(i,j+1);
                    uv_ghost = u0(i,j) * 2 - u0(i+1,j);
                    vu_ghost = v0(i,j) * 2 - v0(i,j+1);
                    vv_ghost = v0(i,j) * 2 - v0(i+1,j);
                    u_wave(i,j) = (aEu(i,j) * u0(i,j+1) + aWu(i,j) * uu_ghost + aSu(i,j)...
                        *uv_ghost + aNu(i,j) * u0(i+1,j) + Suu(i,j)) / aPu(i,j);
                    v_wave(i,j) = (aEv(i,j) * v0(i,j+1) + aWv(i,j) * vu_ghost + aSv(i,j)...
                        *vv_ghost + aNv(i,j) * v0(i+1,j) + Suv(i,j)) / aPv(i,j);
                elseif j > 1 && j < N
                    uv_ghost = u0(i,j) * 2 - u0(i+1,j);
                    vv_ghost = v0(i,j) * 2 - v0(i+1,j);
                    u_wave(i,j) = (aEu(i,j) * u0(i,j+1) + aWu(i,j) * u0(i,j-1) + aSu(i,j)...
                        *uv_ghost + aNu(i,j) * u0(i+1,j) + Suu(i,j)) / aPu(i,j);
                    v_wave(i,j) = (aEv(i,j) * v0(i,j+1) + aWv(i,j) * v0(i,j-1) + aSv(i,j)...
                        *vv_ghost + aNv(i,j) * v0(i+1,j) + Suv(i,j)) / aPv(i,j);
                elseif j == N
                    uu_ghost = u0(i,j) * 2 - u0(i,j-1);
                    uv_ghost = u0(i,j) * 2 - u0(i+1,j);
                    vu_ghost = v0(i,j) * 2 - v0(i,j-1);
                    vv_ghost = v0(i,j) * 2 - v0(i+1,j);
                    u_wave(i,j) = (aEu(i,j) * uu_ghost + aWu(i,j) * u0(i,j-1) + aSu(i,j)...
                        *uv_ghost + aNu(i,j) * u0(i+1,j) + Suu(i,j)) / aPu(i,j);
                    v_wave(i,j) = (aEv(i,j) * vu_ghost + aWv(i,j) * v0(i,j-1) + aSv(i,j)...
                        *vv_ghost + aNv(i,j) * v0(i+1,j) + Suv(i,j)) / aPv(i,j);
                end
            elseif i>1 && i<M
                if j == 1
                    uu_ghost = u0(i,j) * 2 - u0(i,j+1);
                    vu_ghost = v0(i,j) * 2 - v0(i,j+1);
                    u_wave(i,j) = (aEu(i,j) * u0(i,j+1) + aWu(i,j) * uu_ghost + aSu(i,j)...
                        *u0(i-1,j) + aNu(i,j) * u0(i+1,j) + Suu(i,j)) / aPu(i,j);
                    v_wave(i,j) = (aEv(i,j) * v0(i,j+1) + aWv(i,j) * vu_ghost + aSv(i,j)...
                        *v0(i-1,j) + aNv(i,j) * v0(i+1,j) + Suv(i,j)) / aPv(i,j);
                elseif j > 1 && j < N
                    u_wave(i,j) = (aEu(i,j) * u0(i,j+1) + aWu(i,j) * u0(i,j-1) + aSu(i,j)...
                        *u0(i-1,j) + aNu(i,j) * u0(i+1,j) + Suu(i,j)) / aPu(i,j);
                    v_wave(i,j) = (aEv(i,j) * v0(i,j+1) + aWv(i,j) * v0(i,j-1) + aSv(i,j)...
                        *v0(i-1,j) + aNv(i,j) * v0(i+1,j) + Suv(i,j)) / aPv(i,j);
                elseif j == N
                    uu_ghost = u0(i,j) * 2 - u0(i,j-1);
                    vu_ghost = v0(i,j) * 2 - v0(i,j-1);
                    u_wave(i,j) = (aEu(i,j) * uu_ghost + aWu(i,j) * u0(i,j-1) + aSu(i,j)...
                        *u0(i-1,j) + aNu(i,j) * u0(i+1,j) + Suu(i,j)) / aPu(i,j);
                    v_wave(i,j) = (aEv(i,j) * vu_ghost + aWv(i,j) * v0(i,j-1) + aSv(i,j)...
                        *v0(i-1,j) + aNv(i,j) * v0(i+1,j) + Suv(i,j)) / aPv(i,j);
                end
            elseif i == M
                if j == 1
                    uu_ghost = u0(i,j) * 2 - u0(i,j+1);
                    uv_ghost = u0(i,j) * 2 - u0(i-1,j);
                    vu_ghost = v0(i,j) * 2 - v0(i,j+1);
                    vv_ghost = v0(i,j) * 2 - v0(i-1,j);
                    u_wave(i,j) = (aEu(i,j) * u0(i,j+1) + aWu(i,j) * uu_ghost + aSu(i,j)...
                        *u0(i-1,j) + aNu(i,j) * uv_ghost + Suu(i,j)) / aPu(i,j);
                    v_wave(i,j) = (aEv(i,j) * v0(i,j+1) + aWv(i,j) * vu_ghost + aSv(i,j)...
                        *v0(i-1,j) + aNv(i,j) * vv_ghost + Suv(i,j)) / aPv(i,j);
                elseif j > 1 && j < N
                    uv_ghost = u0(i,j) * 2 - u0(i-1,j);
                    vv_ghost = v0(i,j) * 2 - v0(i-1,j);
                    u_wave(i,j) = (aEu(i,j) * u0(i,j+1) + aWu(i,j) * u0(i,j-1) + aSu(i,j)...
                        *u0(i-1,j) + aNu(i,j) * uv_ghost + Suu(i,j)) / aPu(i,j);
                    v_wave(i,j) = (aEv(i,j) * v0(i,j+1) + aWv(i,j) * v0(i,j-1) + aSv(i,j)...
                        *v0(i-1,j) + aNv(i,j) * vv_ghost + Suv(i,j)) / aPv(i,j);
                elseif j == N
                    uu_ghost = u0(i,j) * 2 - u0(i,j-1);
                    uv_ghost = u0(i,j) * 2 - u0(i-1,j);
                    vu_ghost = v0(i,j) * 2 - v0(i,j-1);
                    vv_ghost = v0(i,j) * 2 - v0(i-1,j);
                    u_wave(i,j) = (aEu(i,j) * uu_ghost + aWu(i,j) * u0(i,j-1) + aSu(i,j)...
                        *u0(i-1,j) + aNu(i,j) * uv_ghost + Suu(i,j)) / aPu(i,j);
                    v_wave(i,j) = (aEv(i,j) * vu_ghost + aWv(i,j) * v0(i,j-1) + aSv(i,j)...
                        *v0(i-1,j) + aNv(i,j) * vv_ghost + Suv(i,j)) / aPv(i,j);
                end
            end
%              u_wave(1,1) = (aEu(1,1) * u0(1,2) + aWu(1,1) * (u0(1,1)*2-u0(1,2)) ...
%                  + aSu(1,1) * (u0(1,1)*2-u0(2,1)) + aNu(1,1) * u0(2,1) ...
%                  + Suu(1,1)) / aPu(1,1);
%              v_wave(1,1) = (aEv(1,1) * v0(1,2) + aWv(1,1) * (v0(1,1)*2-v0(1,2)) ...
%                  + aSv(1,1) * (v0(1,1)*2-v0(2,1)) + aNv(1,1) * v0(2,1) ...
%                  + Suv(1,1)) / aPv(1,1);
%              u_wave(1,N) = (aEu(1,N) * (u0(1,N)*2-u0(1,N-1)) + aWu(1,N) * u0(1,N-1)...
%                  + aSu(1,N) *(u0(1,N)*2-u0(2,N)) + aNu(1,N) * u0(2,N) ...
%                  + Suu(1,N)) / aPu(1,N);
%              v_wave(1,N) = (aEv(1,N) * (v0(1,N)*2-v0(1,N-1)) + aWv(1,N) * v0(1,N-1)...
%                  + aSv(1,N) *(v0(1,N)*2-v0(2,N)) + aNv(1,N) * v0(2,N) ...
%                  + Suv(1,N)) / aPv(1,N);
%              u_wave(M,1) = (aEu(M,1) * u0(M,2) + aWu(M,1) * (u0(M,1)*2-u0(M,2)) ...
%                  + aSu(M,1) * u0(M-1,1) + aNu(M,1) * (u0(M,1)*2-u0(M-1,1)) ...
%                  + Suu(M,1)) / aPu(M,1);
%              v_wave(M,1) = (aEv(M,1) * v0(M,2) + aWv(M,1) * (v0(M,1)*2-v0(M,2)) ...
%                  + aSv(M,1) * v0(M-1,1) + aNv(M,1) * (v0(M,1)*2-v0(M-1,1)) ...
%                  + Suv(M,1)) / aPv(M,1);
%              u_wave(M,N) = (aEu(M,N) * (u0(M,N)*2-u0(M,N-1)) + aWu(M,N) * u0(M,N-1) ...
%                  + aSu(M,N) * u0(M-1,1) + aNu(M,N) * (u0(M,N)*2-u0(M-1,N)) ...
%                  + Suu(M,N)) / aPu(M,N);
%              v_wave(M,N) = (aEv(M,N) * (v0(M,N)*2-v0(M,N-1)) + aWv(M,N) * v0(M,N-1) ...
%                  + aSv(M,N) * v0(M-1,1) + aNv(M,N) * (v0(M,N)*2-v0(M-1,N)) ...
%                  + Suv(M,N)) / aPv(M,N);
        end
    end

    %计算d_P
    
    d_P_x = zeros(M,N);
    d_P_y = zeros(M,N);
    for i = 1:M
        for j = 1:N
            d_P_x(i,j) = dy / aPu(i,j);
            d_P_y(i,j) = dx / aPv(i,j);
        end
    end


    %计算界面流速
    for i = 1:M
        for j = 1:N
            if i == 1
                ds_ins = (d_P_y(i,j) * 3 - d_P_y(i+1,j)) / 2;
                dn_ins = (d_P_y(i,j)     + d_P_y(i+1,j)) / 2;
                if j == 1
                    de_ins = (d_P_x(i,j)     + d_P_x(i,j+1)) / 2;
                    dw_ins = (d_P_x(i,j) * 3 - d_P_x(i,j+1)) / 2;

                    vn_wave_ins = (v_wave(i,j) + v_wave(i+1,j)) / 2;
                    vs_wave_ins = (v_wave(i,j) * 3 - v_wave(i+1,j)) / 2;
                    ue_wave_ins = (u_wave(i,j) + u_wave(i,j+1)) / 2;
                    uw_wave_ins = (u_wave(i,j) * 3 - u_wave(i,j+1)) / 2;

                    vn(i,j) = vn_wave_ins - dn_ins...
                        * (p0(i+1,j) - p0(i,j));
                    vs(i,j) = 0;
                    ue(i,j) = ue_wave_ins - de_ins...
                        * (p0(i,j+1) - p0(i,j));
                    uw(i,j) = 0;
                elseif j > 1 && j < N
                    de_ins = (d_P_x(i,j)     + d_P_x(i,j+1)) / 2;
                    dw_ins = (d_P_x(i,j)     + d_P_x(i,j-1)) / 2;

                    vn_wave_ins = (v_wave(i,j) + v_wave(i+1,j)) / 2;
                    vs_wave_ins = 0;
                    ue_wave_ins = (u_wave(i,j) + u_wave(i,j+1)) / 2;
                    uw_wave_ins = (u_wave(i,j) + u_wave(i,j-1)) / 2;

                    vn(i,j) = vn_wave_ins - dn_ins...
                        * (p0(i+1,j) - p0(i,j));
                    vs(i,j) = 0;
                    ue(i,j) = ue_wave_ins - de_ins...
                        * (p0(i,j+1) - p0(i,j));
                    uw(i,j) = uw_wave_ins - dw_ins...
                        * (p0(i,j) - p0(i,j-1));
                elseif j == N
                    de_ins = (d_P_x(i,j) * 3 - d_P_x(i,j-1)) / 2;
                    dw_ins = (d_P_x(i,j)     + d_P_x(i,j-1)) / 2;

                    vn_wave_ins = (v_wave(i,j) + v_wave(i+1,j)) / 2;
                    vs_wave_ins = (v_wave(i,j) * 3 - v_wave(i+1,j)) / 2;
                    ue_wave_ins = (u_wave(i,j) * 3 - u_wave(i,j-1)) / 2;
                    uw_wave_ins = (u_wave(i,j) + u_wave(i,j-1)) / 2;

                    vn(i,j) = vn_wave_ins - dn_ins * ...
                        (p0(i+1,j) - p0(i,j));
                    vs(i,j) = 0;
                    ue(i,j) = 0;
                    uw(i,j) = uw_wave_ins - dw_ins *...
                        (p0(i,j) - p0(i,j-1));
                end
            elseif i > 1 && i < M
                if j == 1
                    dn_ins = (d_P_y(i,j)     + d_P_y(i+1,j)) / 2;
                    ds_ins = (d_P_y(i,j)     + d_P_y(i-1,j)) / 2;
                    de_ins = (d_P_x(i,j)     + d_P_x(i,j+1)) / 2;
                    dw_ins = 0;

                    vn_wave_ins = (v_wave(i,j) + v_wave(i+1,j)) / 2;
                    vs_wave_ins = (v_wave(i,j) + v_wave(i-1,j)) / 2;
                    ue_wave_ins = (u_wave(i,j) + u_wave(i,j+1)) / 2;
                    uw_wave_ins = 0;

                    vn(i,j) = vn_wave_ins - dn_ins * ...
                        (p0(i+1,j) - p0(i,j));
                    vs(i,j) = vs_wave_ins - ds_ins *...
                        (p0(i,j) - p0(i-1,j));
                    ue(i,j) = ue_wave_ins - de_ins *...
                        (p0(i,j+1) - p0(i,j));
                    uw(i,j) = 0;
                elseif j > 1 && j < N
                    dn_ins = (d_P_y(i,j)     + d_P_y(i+1,j)) / 2;
                    ds_ins = (d_P_y(i,j)     + d_P_y(i-1,j)) / 2;
                    de_ins = (d_P_x(i,j)     + d_P_x(i,j+1)) / 2;
                    dw_ins = (d_P_x(i,j)     + d_P_x(i,j-1)) / 2;

                    vn_wave_ins = (v_wave(i,j) + v_wave(i+1,j)) / 2;
                    vs_wave_ins = (v_wave(i,j) + v_wave(i-1,j)) / 2;
                    ue_wave_ins = (u_wave(i,j) + u_wave(i,j+1)) / 2;
                    uw_wave_ins = (u_wave(i,j) + u_wave(i,j-1)) / 2;

                    vn(i,j) = vn_wave_ins - dn_ins * ...
                        (p0(i+1,j) - p0(i,j));
                    vs(i,j) = vs_wave_ins - ds_ins *...
                        (p0(i,j) - p0(i-1,j));
                    ue(i,j) = ue_wave_ins - de_ins *...
                        (p0(i,j+1) - p0(i,j));
                    uw(i,j) = uw_wave_ins - dw_ins *...
                        (p0(i,j) - p0(i,j-1));
                elseif j == N
                    dn_ins = (d_P_y(i,j)     + d_P_y(i+1,j)) / 2;
                    ds_ins = (d_P_y(i,j)     + d_P_y(i-1,j)) / 2;
                    de_ins = 0;
                    dw_ins = (d_P_x(i,j)     + d_P_x(i,j-1)) / 2;

                    vn_wave_ins = (v_wave(i,j) + v_wave(i+1,j)) / 2;
                    vs_wave_ins = (v_wave(i,j) + v_wave(i-1,j)) / 2;
                    ue_wave_ins = 0;
                    uw_wave_ins = (u_wave(i,j) + u_wave(i,j-1)) / 2;

                    vn(i,j) = vn_wave_ins - dn_ins * ...
                        (p0(i+1,j) - p0(i,j));
                    vs(i,j) = vs_wave_ins - ds_ins *...
                        (p0(i,j) - p0(i-1,j));
                    ue(i,j) = 0;
                    uw(i,j) = uw_wave_ins - dw_ins *...
                        (p0(i,j) - p0(i,j-1));
                end
            elseif i == M
                if j == 1
                    dn_ins = 0;
                    ds_ins = (d_P_y(i,j)     + d_P_y(i-1,j)) / 2;
                    de_ins = (d_P_x(i,j)     + d_P_x(i,j+1)) / 2;
                    dw_ins = 0;

                    vn_wave_ins = 0;
                    vs_wave_ins = (v_wave(i,j) + v_wave(i-1,j)) / 2;
                    ue_wave_ins = (u_wave(i,j) + u_wave(i,j+1)) / 2;
                    uw_wave_ins = 0;

                    vn(i,j) = 0;
                    vs(i,j) = vs_wave_ins - ds_ins *...
                        (p0(i,j) - p0(i-1,j));
                    ue(i,j) = ue_wave_ins - de_ins *...
                        (p0(i,j+1) - p0(i,j));
                    uw(i,j) = 0;
                elseif j > 1 && j < N
                    dn_ins = 0;
                    ds_ins = (d_P_y(i,j)     + d_P_y(i-1,j)) / 2;
                    de_ins = (d_P_x(i,j)     + d_P_x(i,j+1)) / 2;
                    dw_ins = (d_P_x(i,j)     + d_P_x(i,j-1)) / 2;

                    vn_wave_ins = 0;
                    vs_wave_ins = (v_wave(i,j) + v_wave(i-1,j)) / 2;
                    ue_wave_ins = (u_wave(i,j) + u_wave(i,j+1)) / 2;
                    uw_wave_ins = (u_wave(i,j) + u_wave(i,j-1)) / 2;

                    vn(i,j) = 0;
                    vs(i,j) = vs_wave_ins - ds_ins *...
                        (p0(i,j) - p0(i-1,j));
                    ue(i,j) = ue_wave_ins - de_ins *...
                        (p0(i,j+1) - p0(i,j));
                    uw(i,j) = uw_wave_ins - dw_ins *...
                        (p0(i,j) - p0(i,j-1));
                elseif j == N
                    dn_ins = 0;
                    ds_ins = (d_P_y(i,j)     + d_P_y(i-1,j)) / 2;
                    de_ins = 0;
                    dw_ins = (d_P_x(i,j)     + d_P_x(i,j-1)) / 2;

                    vn_wave_ins = 0;
                    vs_wave_ins = (v_wave(i,j) + v_wave(i-1,j)) / 2;
                    ue_wave_ins = 0;
                    uw_wave_ins = (u_wave(i,j) + u_wave(i,j-1)) / 2;

                    vn(i,j) = 0;
                    vs(i,j) = vs_wave_ins - ds_ins *...
                        (p0(i,j) - p0(i-1,j));
                    ue(i,j) = 0;
                    uw(i,j) = uw_wave_ins - dw_ins *...
                        (p0(i,j) - p0(i,j-1));
                end
            end
        end
    end
    

    %calculate p_cor's source:b
    %     for j = 1:N
    %         for i = 1:M
    %             if i>1 && i<M && j<N && j>1
    %                 b(i,j) = (rho(i,j-1)+rho(i,j)) /2 * uw(i,j) * dy ...
    %                     -(rho(i,j)+rho(i,j+1)) /2 * ue(i,j) * dy ...
    %                     +(rho(i-1,j)+rho(i,j)) /2 * vs(i,j) * dx ...
    %                     -(rho(i,j)+rho(i+1,j)) /2 * vn(i,j) * dx;
    %             end
    %             if i == 1 && j ~= 1 && j ~= N
    %                 b(i,j) = (rho(i,j-1)+rho(i,j)) /2 * uw(i,j) * dy ...
    %                     -(rho(i,j)+rho(i,j+1)) /2 * ue(i,j) * dy ...
    %                     +(rho(i,j)*3-rho(i+1,j)) /2 * vs(i,j) * dx ...
    %                     -(rho(i,j)+rho(i+1,j)) /2 * vn(i,j) * dx;
    %             end
    %             if i == M && j ~= 1 && j ~= N
    %                 b(i,j) = (rho(i,j-1)+rho(i,j)) /2 * uw(i,j) * dy ...
    %                     -(rho(i,j)+rho(i,j+1)) /2 * ue(i,j) * dy ...
    %                     +(rho(i-1,j)+rho(i,j)) /2 * vs(i,j) * dx ...
    %                     -(rho(i,j)*3-rho(i-1,j)) /2 * vn(i,j) * dx;
    %             end
    %             if j == 1 && i ~= 1 && i ~= M
    %                 b(i,j) = (rho(i,j)*3-rho(i,j+1)) /2 * uw(i,j) * dy ...
    %                     -(rho(i,j)+rho(i,j+1)) /2 * ue(i,j) * dy ...
    %                     +(rho(i-1,j)+rho(i,j)) /2 * vs(i,j) * dx ...
    %                     -(rho(i,j)+rho(i+1,j)) /2 * vn(i,j) * dx;
    %             end
    %             if j == M && i ~= 1 && i ~= M
    %                 b(i,j) = (rho(i,j-1)+rho(i,j)) /2 * uw(i,j) * dy ...
    %                     -(rho(i,j)*3-rho(i,j-1)) /2 * ue(i,j) * dy ...
    %                     +(rho(i-1,j)+rho(i,j)) /2 * vs(i,j) * dx ...
    %                     -(rho(i,j)+rho(i+1,j)) /2 * vn(i,j) * dx;
    %             end
    %         end
    %     end
    %     for i = 1:M
    %         for j = 1:N
    %             if i >1
    %                 ds(i,j) = dy * 0.5 * (1/aPv(i,j) + 1 / aPv(i-1,j));
    %             end
    %             if i == 1
    %                 ds(i,j) = dy * 0.5 * ( 1/aPv(i,j) + 1/(aPv(i,j) * 2 - aPv(i+1,j)));
    %             end
    %             if i < M
    %                 dn(i,j) = dy * 0.5 * (1/aPv(i,j) + 1 / aPv(i+1,j));
    %             end
    %             if i == M
    %                 dn(i,j) = dy * 0.5 * ( 1/aPv(i,j) + 1/(aPv(i,j) * 2 - aPv(i-1,j)));
    %             end
    %             if j < M
    %                 de(i,j) = dx * 0.5 * (1/aPu(i,j) + 1 / aPu(i,j+1));
    %             end
    %             if j == M
    %                 de(i,j) = dx * 0.5 * (1/aPu(i,j) + 1/(aPu(i,j) * 2 - aPu(i,j-1)));
    %             end
    %             if j > 1
    %                 dw(i,j) = dx * 0.5 * (1/aPu(i,j) + 1 / aPu(i,j-1));
    %             end
    %             if j == 1
    %                 dw(i,j) = dx * 0.5 * (1/aPu(i,j) + 1 / (aPu(i,j) * 2 - aPu(i,j+1)));
    %             end
    %         end
    %     end
    %
    %     aE = rho * de ;
    %     aW = rho * dw ;
    %     aN = rho * dn ;
    %     aS = rho * ds ;
    %     aP = aE + aN + aS + aW;
    %     aE(1:M,N) = zeros(M,1);
    %     aN(1,1:N) = zeros(1,N);
    %     aS(M,1:N) = zeros(1,N);
    %     aW(1:M,1) = zeros(M,1);
    ppp = 0;
    r = zeros(M*N + M*(N-1)*2 + N*(M-1)*2,1);
    c = zeros(M*N + M*(N-1)*2 + N*(M-1)*2,1);
    e = zeros(M*N + M*(N-1)*2 + N*(M-1)*2,1);
    for j = 1:N
        for i = 1:M
            uw_ins = uw(i,j);
            ue_ins = ue(i,j);
            vn_ins = vn(i,j);
            vs_ins = vs(i,j);
            if i == 1
                rho_s_ins = (rho(i,j) * 3 - rho(i+1,j)) / 2;
                rho_n_ins = (rho(i+1,j) + rho(i,j)) / 2;

                aP_n_ins  = (1/aPv(i,j) + 1/aPv(i+1,j)) / 2;
                aP_s_ins  = (1/aPv(i,j) + 1 / (aPv(i,j) * 2- aPv(i+1,j))) / 2;
                if j == 1
                    rho_w_ins = (rho(i,j) * 3 - rho(i,j+1)) / 2;
                    rho_e_ins = (rho(i,j+1) + rho(i,j)) / 2;

                    aP_w_ins  = (1/aPu(i,j) + 1 / (aPu(i,j) * 2 - aPu(i,j+1))) / 2;
                    aP_e_ins  = (1/aPu(i,j) + 1 / aPu(i,j+1)) / 2;
                elseif j>1 && j<N
                    rho_w_ins = (rho(i,j-1) + rho(i,j)) / 2;
                    rho_e_ins = (rho(i,j+1) + rho(i,j)) / 2;

                    aP_w_ins  = (1/aPu(i,j) + 1 / aPu(i,j-1)) / 2;
                    aP_e_ins  = (1/aPu(i,j) + 1 / aPu(i,j+1)) / 2;
                elseif j == N
                    rho_w_ins = (rho(i,j-1) + rho(i,j)) / 2;
                    rho_e_ins = (rho(i,j) * 3 - rho(i,j-1)) / 2;

                    aP_w_ins  = (1/aPu(i,j) + 1 / aPu(i,j-1)) / 2;
                    aP_e_ins  = (1/aPu(i,j) + 1/(aPu(i,j) * 2 - aPu(i,j-1))) / 2;
                end
            elseif i>1 && i<M
                rho_s_ins = (rho(i-1,j) + rho(i,j)) / 2;
                rho_n_ins = (rho(i+1,j) + rho(i,j)) / 2;

                aP_n_ins  = (1/aPv(i,j) + 1 / aPv(i+1,j)) / 2;
                aP_s_ins  = (1/aPv(i,j) + 1 / aPv(i-1,j)) / 2;
                if j == 1
                    rho_w_ins = (rho(i,j) * 3 - rho(i,j+1)) / 2;
                    rho_e_ins = (rho(i,j+1) + rho(i,j)) / 2;

                    aP_w_ins  = (1/aPu(i,j) + 1 / (aPu(i,j) * 2 - aPu(i,j+1))) / 2;
                    aP_e_ins  = (1/aPu(i,j) + 1 / aPu(i,j+1)) / 2;
                elseif j>1 && j<N
                    rho_w_ins = (rho(i,j-1) + rho(i,j)) / 2;
                    rho_e_ins = (rho(i,j+1) + rho(i,j)) / 2;

                    aP_w_ins  = (1/aPu(i,j) + 1 / aPu(i,j-1)) / 2;
                    aP_e_ins  = (1/aPu(i,j) + 1 / aPu(i,j+1)) / 2;
                elseif j == N
                    rho_w_ins = (rho(i,j-1) + rho(i,j)) / 2;
                    rho_e_ins = (rho(i,j) * 3 - rho(i,j-1)) / 2;

                    aP_w_ins  = (1/aPu(i,j) + 1 / aPu(i,j-1)) / 2;
                    aP_e_ins  = (1/aPu(i,j) + 1/(aPu(i,j) * 2 - aPu(i,j-1))) / 2;
                end
            elseif i == M
                rho_s_ins = (rho(i-1,j) + rho(i,j)) / 2;
                rho_n_ins = (rho(i,j) * 3 - rho(i-1,j)) / 2;

                aP_n_ins  = (1/aPv(i,j) + 1 / (aPv(i,j) * 2 - aPv(i-1,j))) / 2;
                aP_s_ins  = (1/aPv(i,j) + 1 / aPv(i-1,j)) / 2;
                if j == 1
                    rho_w_ins = (rho(i,j) * 3 - rho(i,j+1)) / 2;
                    rho_e_ins = (rho(i,j+1) + rho(i,j)) / 2;

                    aP_w_ins  = (1/aPu(i,j) + 1 / (aPu(i,j) * 2 - aPu(i,j+1))) / 2;
                    aP_e_ins  = (1/aPu(i,j) + 1 / aPu(i,j+1)) / 2;
                elseif j>1 && j<N
                    rho_w_ins = (rho(i,j-1) + rho(i,j)) / 2;
                    rho_e_ins = (rho(i,j+1) + rho(i,j)) / 2;

                    aP_w_ins  = (1/aPu(i,j) + 1 / aPu(i,j-1)) / 2;
                    aP_e_ins  = (1/aPu(i,j) + 1 / aPu(i,j+1)) / 2;
                elseif j == N
                    rho_w_ins = (rho(i,j-1) + rho(i,j)) / 2;
                    rho_e_ins = (rho(i,j) * 3 - rho(i,j-1)) / 2;

                    aP_w_ins  = (1/aPu(i,j) + 1 / aPu(i,j-1)) / 2;
                    aP_e_ins  = (1/aPu(i,j) + 1/(aPu(i,j) * 2 - aPu(i,j-1))) / 2;
                end
            end
            b(i,j) = rho_w_ins * uw_ins * dy - rho_e_ins * ue_ins * dy +...
                rho_s_ins * vs_ins * dx - rho_n_ins * vn_ins * dx;

            aP_W = rho_w_ins * dx * dx / aP_w_ins;
            aP_E = rho_e_ins * dx * dx / aP_e_ins;
            aP_S = rho_s_ins * dy * dy / aP_s_ins;
            aP_N = rho_n_ins * dy * dy / aP_n_ins;
%             if i == 1
%                 aP_S = 1e10;
%             end
%             if i == 1 && j == N
%                 aP_S = 1e10;
%             end
%             if j == 1
%                 aP_W = 1e10;
%             end
%             if j == N
%                 aP_E = 1e10;
%             end
            aP_P = aP_W + aP_E + aP_S + aP_N;
            I = (j - 1) * M + i;
            if I
                ppp = ppp + 1;
                e(ppp) = aP_P;
                r(ppp) = I ;
                c(ppp) = I ;
            end
            if I<M*N
                ppp = ppp + 1;
                e(ppp) = -aP_S;
                r(ppp) = I ;
                c(ppp) = I+1 ;
            end
            if I<=M*(N-1)
                ppp = ppp + 1;
                e(ppp) = -aP_E;
                r(ppp) = I ;
                c(ppp) = I+M ;
            end
            if I>M
                ppp = ppp + 1;
                e(ppp) = -aP_W;
                r(ppp) = I ;
                c(ppp) = I-M ;
            end
            if I>1
                ppp = ppp + 1;
                e(ppp) = -aP_N;
                r(ppp) = I ;
                c(ppp) = I-1 ;
            end
        end
    end

    matrix = sparse(r,c,e);
    b = reshape_a(b,M,N);

    back_p_cor = zeros(M,N);
    p_cor = matrix \ b;
    for i = 1 : M
        for j=1:N
            back_p_cor(i,j) = p_cor((j-1)*M+i,1);%dimension up
        end
    end

    back_b = zeros(M,N);
    
    for j = 1 : N
        for i = 1 : M
            back_b(i,j) = b((j-1)*M+i,1);%dimension up
        end
    end
    
    p_cor = back_p_cor;
    b = back_b;

    %calculat u_cor & v_cor
    for i = 1:M
        for j = 1:N
            if j < N && j > 1
                u_cor(i,j) = dx / aPu(i,j) * (p_cor(i,j-1) - p_cor(i,j+1))/2;
            end
            if j == N
                p_ghost = p_cor(i,j) * 2 - p_cor(i,j-1);
                u_cor(i,j) = dx / aPu(i,j) * (p_cor(i,j-1) - p_ghost)/2;
            end
            if j == 1
                p_ghost = p_cor(i,j) * 2 - p_cor(i,j+1);
                u_cor(i,j) = dx / aPu(i,j) * (p_ghost - p_cor(i,j+1))/2;
            end
            if i > 1 && i<M
                v_cor(i,j) = dy / aPv(i,j) * (p_cor(i-1,j) - p_cor(i+1,j))/2;
            end
            if i == 1
                p_ghost = p_cor(i,j) * 2 - p_cor(i+1,j);
                v_cor(i,j) = dy / aPv(i,j) * (p_ghost - p_cor(i+1,j))/2;
            end
            if i == M
                p_ghost = p_cor(i,j) * 2 - p_cor(i-1,j);
                v_cor(i,j) = dy / aPv(i,j) * (p_cor(i-1,j) - p_ghost)/2;
            end
        end
    end

%     for j = 1:N
%         for i = 1:M
%             if i == M
%                 u(i,j) = uN(y(j));
%             end
%             if i == 1
%                 u(i,j) = uS(y(j));
%             end
%             if j == 1
%                 v(i,j) = vW(x(i));
%             end
%             if j == N
%                 v(i,j) = vE(x(i));
%             end
%         end
%     end

    p_cor(1,N) = 0;
    p = p0 + alpha_p * p_cor;
    u = u0 + u_cor;
    v = v0 + v_cor;
    n = n + 1;
    norm_b(n)  = norm(b,2);
    norm_v0(n) = norm(v_cor,2);
    norm_u0(n) = norm(u_cor,2);
    k = n;
    if mod(k,1000) == 0
%                 figure(1);
%                 subplot(4,4,1);
%                 mesh(y,x,u);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('u');
%                 box on;
%                 title('u');
%                 subplot(4,4,2);
%                 mesh(y,x,v);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('v');
%                 box on;
%                 title('v');
%                 subplot(4,4,3);
%                 mesh(y,x,p0);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('p');
%                 box on;
%                 title('p0');
%                 subplot(4,4,5);
%                 mesh(y,x,u_wave);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('u wave');
%                 box on;
%                 title('u wave');
%                 subplot(4,4,6);
%                 mesh(y,x,v_wave);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('v wave');
%                 box on;
%                 title('v wave');
%                 subplot(4,4,9);
%                 mesh(y,x,ue);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('ue');
%                 box on;
%                 title('ue');
%                 subplot(4,4,4);
%                 matrix = full(matrix);
%                 dxx = 1/M/N;
%                 xx  = [0:M*N-1] * dxx;
%                 mesh(xx,xx,matrix);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('mtrx');
%                 box on;
%                 title('matrix');
%                 subplot(4,4,10);
%                 mesh(x,y,uw);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('uw');
%                 box on;
%                 title('uw');
%                 subplot(4,4,11);
%                 mesh(x,y,vn);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('vn');
%                 title('vn');
%                 subplot(4,4,12);
%                 mesh(x,y,vs);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('vn');
%                 title('vs');
%                 subplot(4,4,13);
%                 mesh(x,y,u_cor);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('u cor');
%                 title('u cor');
%                 subplot(4,4,14);
%                 mesh(x,y,v_cor);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('v cor');
%                 title('v cor');
%                 subplot(4,4,15);
%                 mesh(x,y,p_cor);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('p cor');
%                 title('p cor');
%                 subplot(4,4,16);
%                 mesh(x,y,b);
%                 xlabel('x(m)');
%                 ylabel('y(m)');
%                 zlabel('b');
%                 title('b');
%                 subplot(4,4,7);
%                 mesh(y,x,aPu);
%                 xlabel('x');
%                 ylabel('y');
%                 title('aPu');
%                 subplot(4,4,8);
%                 mesh(y,x,aPv);
%                 xlabel('x');
%                 ylabel('y');
%                 title('aPv');
                figure(2);
                plot(1:n,norm_u0(1:n),'g',1:n,norm_v0(1:n),'r',1:n,norm_b(1:n),'b');
                %     semilogy(n,10.^(-n));
                legend('u0','v0','b');

    
    if norm_b(n)<= 1e-6 && n >= 10
        break;
    end
    

    if mod(k,10000) == 0
    clc;
    fprintf('\niteration = %d',n);
    fprintf('\nnorm b = %d',norm_b(n));
    end
    end
end
end
