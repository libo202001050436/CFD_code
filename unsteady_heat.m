function [T,u,v] = unsteady_heat(Lx,Ly,time,dt,T0,TW,TE,TN,TS,dx,dy,rho0,...
    alpha,V_char,L_char,T_char,mu,g,P_char,heat_cap,k,uN,uS,vW,vE)

theta = 0.5;
M1 = Lx / dx;
x  = dx / 2 + [0:M1-1] * dx;

M2 = Ly / dy;
y  = dy / 2 + [0:M2-1] * dy;

Mt = time / dt;

u = zeros(M1,M2);
v = zeros(M1,M2);
T = zeros(M1,M2);

for j = 1:M2
    for i = 1:M1
        T(i,j) = T0;
    end
end
Sp = zeros(M1,M2);
Su = zeros(M1,M2);
rho= zeros(M1,M2);
Eu = zeros(M1,M2);
Re = zeros(M1,M2);


for j = 1:M2
    for i = 1:M1
        rho(i,j) = rho0;
    end
end
% %%
% rho(2*M1/5:3*M1/5,2*M2/5:3*M2/5) = 0.9;
% %%
%%在NS的基础上加入能量方程
for t = 0:Mt
    %计算温度场
    if t>0
        u = u0;
        v = v0;
    end
    rho0 = rho;
    ppp=0;
    %%
%     for i = 1:M1
%         for j = 1:M2
%             kapa = k / rho(i,j) / heat_cap;
%             coef_ins = kapa * T_char * theta * dt / L_char / L_char;
%             u_ins = u(i,j);
%             v_ins = v(i,j);
%             TP_ins  = T(i,j);
%             aP0 = (dx * dy / dt - kapa * T_char / L_char / L_char...
%                 * (1 - theta) * 2 * dy * dt / dx -...
%                 kapa * T_char / L_char / L_char * (1 - theta) ...
%                 * 2 * dx * dt / dy);
%             aE0 = (-(1 - theta) * v_ins * dx * dt / 2 + kapa * T_char...
%                 /L_char / L_char * (1-theta) * dy * dt / dy);
%             aW0 = aE0;
%             aN0 = (-(1 - theta) * v_ins * dx * dt / 2 + kapa * T_char...
%                 /L_char / L_char * (1-theta) * dy * dt / dy);
%             aS0 = aN0;
%             if i>1 && i<M1
%                 TN_ins = T(i+1,j);
%                 TS_ins = T(i-1,j);
%                 if j>1 && j<M2
%                     TE_ins = T(i,j+1);
%                     TW_ins = T(i,j-1);
%                     aW(i,j) = theta * (coef_ins * dy / dx + theta * u_ins *dy*dt/2);
%                     aE(i,j) = theta * (coef_ins * dy / dx - theta * u_ins *dy*dt/2);
%                     aS(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     aN(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     Sp(i,j) = -(dx * dy / dt);
%                 end
%                 if j == 1
%                     TW_ins = T0;
%                     TE_ins = T(i,j+1);
%                     aE(i,j) = theta * (coef_ins * dy / dx - theta * u_ins *dy*dt/2);
%                     aS(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     aN(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     Sp(i,j) = -(dx * dy / dt + theta * u_ins * dy * dt ...
%                         + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
%                 end
%                 if j == M2
%                     TE_ins = T0;
%                     TW_ins = T(i,j-1);
%                     aW(i,j) = theta * (coef_ins * dy / dx + theta * u_ins *dy*dt/2);
%                     aS(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     aN(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     Sp(i,j) = -(dx * dy / dt - theta * u_ins * dy * dt ...
%                         + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
%                 end
%             end
%             if i == 1
%                 TN_ins = T(i+1,j);
%                 TS_ins = TS(x(i));
%                 if j<M2 && j>1
%                     TE_ins = T(i,j+1);
%                     TW_ins = T(i,j-1);
%                     aW(i,j) = theta * (coef_ins * dy / dx + theta * u_ins *dy*dt/2);
%                     aE(i,j) = theta * (coef_ins * dy / dx - theta * u_ins *dy*dt/2);
%                     aN(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     Sp(i,j) = -(dx * dy / dt + theta * v_ins * dx * dt ...
%                         + 2 * kapa * T_char * theta *dx*dt / dy / L_char / L_char);
%                 end
%                 if j == 1
%                     TE_ins = T(i,j+1);
%                     TW_ins = T0;
%                     aE(i,j) = theta * (coef_ins * dy / dx - theta * u_ins *dy*dt/2);
%                     aN(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     Sp(i,j) = -(dx * dy / dt + theta * v_ins * dy * dt ...
%                         + 2 * kapa * T_char * theta *dx*dt / ...
%                         dy / L_char / L_char + theta * u_ins * dy * dt ...
%                         + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
%                 end
%                 if j == M2
%                     TE_ins = T0;
%                     TW_ins = T(i,j-1);
%                     aW(i,j) = theta * (coef_ins * dy / dx + theta * u_ins *dy*dt/2);
%                     aN(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     Sp(i,j) = -(dx * dy / dt + theta * v_ins * dy * dt ...
%                         + 2 * kapa * T_char * theta *dy*dt / ...
%                         dx / L_char / L_char - theta * u_ins * dy * dt ...
%                         + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
%                 end
%             end
%             if i == M1
%                 TN_ins = T0;
%                 TS_ins = T(i-1,j);
%                 if j<M2 && j>1
%                     TE_ins = T(i,j+1);
%                     TW_ins = T(i,j-1);
%                     aW(i,j) = theta * (coef_ins * dy / dx + theta * u_ins *dy*dt/2);
%                     aE(i,j) = theta * (coef_ins * dy / dx - theta * u_ins *dy*dt/2);
%                     aS(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     Sp(i,j) = -(dx * dy / dt - theta * v_ins * dx * dt ...
%                         + 2 * kapa * T_char * theta *dx*dt / dy / L_char / L_char);
%                 end
%                 if j == 1
%                     TE_ins = T(i,j+1);
%                     TW_ins = T0;
%                     aE(i,j) = theta * (coef_ins * dy / dx - theta * u_ins *dy*dt/2);
%                     aS(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     Sp(i,j) = -(dx * dy / dt - theta * v_ins * dx * dt ...
%                         + 2 * kapa * T_char * theta *dx*dt / dy / L_char ...
%                         / L_char + theta * u_ins * dy * dt ...
%                         + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
%                 end
%                 if j == M2
%                     TE_ins = T(i,j) * 2 - T(i,j-1);
%                     TW_ins = T(i,j-1);
%                     aW(i,j) = theta * (coef_ins * dy / dx + theta * u_ins *dy*dt/2);
%                     aS(i,j) = theta * (coef_ins * dx / dy + theta * v_ins *dx*dt/2);
%                     Sp(i,j) = -(dx * dy / dt - theta * v_ins * dx * dt ...
%                         + 2 * kapa * T_char * theta *dx*dt / dy / L_char ...
%                         / L_char - theta * u_ins * dy * dt ...
%                         + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
%                 end
%             end
%             Su(i,j) = aP0 * TP_ins + aE0 * TE_ins + aW0 * TW_ins...
%                 + aN0 * TN_ins + aS0 * TS_ins;
%         end
%     end
%     aS(1,1:M2) = zeros(1,M2);
%     aN(M1,1:M2)= zeros(1,M2);
%     aW(1:M1,1) = zeros(M1,1);
%     aE(1:M1,M2)= zeros(M1,1);
%     aP = aW + aN + aS + aE - Sp;
%     T1 = zeros(M1,M2);
%     T2 = T1;
%     matrix = sparse_coef_auto(aP,aW,aE,aN,aS,M1,M2);
%     %     re_Su  = reshape_a(Su,M1,M2);
%     %     T_downD = matrix \ re_Su;
%     %     for j = 1 : M2
%     %         for i = 1 : M1
%     %             T1(i,j) = T_downD((j-1)*M1+i,1);
%     %         end
%     %     end
%     re_Su  = reshape_a(Su,M1,M2);
%     T_downD = matrix \ re_Su;
%     for j = 1 : M2
%         for i = 1 : M1
%             T2(i,j) = T_downD((j-1)*M1+i,1);
%         end
%     end
%     %     T2 = T2;
%     T = T2;
%     for j = 1:M2
%         for i = 1:M1
%             rho(i,j) = (1 - alpha * T(i,j)) * rho0(i,j);
%         end
%     end
    %%
    r = zeros(M1*M2 + M1*(M2-1)*2 + M2*(M1-1)*2,1);
    c = zeros(M1*M2 + M1*(M2-1)*2 + M2*(M1-1)*2,1);
    e = zeros(M1*M2 + M1*(M2-1)*2 + M2*(M1-1)*2,1);
    for j = 1:M2
        for i = 1:M1
            rho_ins = rho0(i,j);
            kapa = k / rho_ins / heat_cap;
            u_ins = u(i,j);
            v_ins = v(i,j);
            T_P_ins = T(i,j);
            aP0 = (dx * dy / dt - kapa * T_char / L_char / L_char...
                * (1 - theta) * 2 * dy * dt / dx -...
                kapa * T_char / L_char / L_char * (1 - theta) ...
                * 2 * dx * dt / dy);
            aE0 = (-(1 - theta) * v_ins * dx * dt / 2 + kapa * T_char...
                /L_char / L_char * (1-theta) * dy * dt / dy);
            aW0 = aE0;
            aN0 = (-(1 - theta) * v_ins * dx * dt / 2 + kapa * T_char...
                /L_char / L_char * (1-theta) * dy * dt / dy);
            aS0 = aN0;
            L2 = L_char * L_char;
            if i == 1
                aN_ins = kapa * T_char * theta * dx * dt / L2 / dy - ...
                    theta * v_ins * dx * dt / 2;
                aS_ins = 0;
                T_N_ins = T(i+1,j);
%                 if t == 0
                    T_S_ins = TS(y(j));
%                 elseif t
%                     T_S_ins = T(i,j) * 2 - T(i+1,j);
%                 end
                if j == 1
                    aW_ins = 0;
                    aE_ins = kapa * T_char * theta * dy * dt / L2 / dx - ...
                        theta * u_ins * dy * dt / 2;
                    T_E_ins = T(i,j+1);
                    T_W_ins = TW(y(i));
                    Sp = -(dx * dy / dt + theta * v_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dx*dt / ...
                        dy / L_char / L_char + theta * u_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
                elseif j<M2 && j>1
                    aW_ins = kapa * T_char * theta * dy * dt / L2 / dx + ...
                        theta * u_ins * dy * dt / 2;
                    aE_ins = kapa * T_char * theta * dy * dt / L2 / dx - ...
                        theta * u_ins * dy * dt / 2;
                    T_E_ins = T(i,j+1);
                    T_W_ins = T(i,j-1);
                    Sp = -(dx * dy / dt + theta * v_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dx*dt / ...
                        dy / L_char / L_char + theta * u_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
                elseif j == M2
                    aE_ins = 0;
                    aW_ins = kapa * T_char * theta * dy * dt / L2 / dx + ...
                        theta * u_ins * dy * dt / 2;
                    T_E_ins = TE(y(i));
                    T_W_ins = T(i,j-1);
                    Sp = -(dx * dy / dt + theta * v_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dy*dt / ...
                        dx / L_char / L_char - theta * u_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
                end
            elseif i>1 && i<M1
                aN_ins = kapa * T_char * theta * dx * dt / L2 / dy - ...
                    theta * v_ins * dx * dt / 2;
                aS_ins = kapa * T_char * theta * dx * dt / L2 / dy + ...
                    theta * v_ins * dx * dt / 2;
                T_N_ins = T(i+1,j);
                T_S_ins = T(i-1,j);
                if j == 1
                    aW_ins = 0;
                    aE_ins = kapa * T_char * theta * dy * dt / L2 / dx - ...
                        theta * u_ins * dy * dt / 2;
                    T_E_ins = T(i,j+1);
                    T_W_ins = TW(y(i));
                    Sp = -(dx * dy / dt + theta * u_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dy*dt / dx / L_char / L_char);
                elseif j<M2 && j>1
                    aW_ins = kapa * T_char * theta * dy * dt / L2 / dx + ...
                        theta * u_ins * dy * dt / 2;
                    aE_ins = kapa * T_char * theta * dy * dt / L2 / dx - ...
                        theta * u_ins * dy * dt / 2;
                    T_E_ins = T(i,j+1);
                    T_W_ins = T(i,j-1);
                    Sp = -(dx * dy / dt);
                elseif j == M2
                    aW_ins = kapa * T_char * theta * dy * dt / L2 / dx + ...
                        theta * u_ins * dy * dt / 2;
                    aE_ins = 0;
                    T_E_ins = TE(y(i));
                    T_W_ins = T(i,j-1);
                    Sp = -(dx * dy / dt - theta * u_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dy*dt / dx ...
                        / L_char / L_char);
                end
            elseif i == M1
                aN_ins = 0;
                aS_ins = kapa * T_char * theta * dx * dt / L2 / dy + ...
                    theta * v_ins * dx * dt / 2;
                Sp = -(dx * dy / dt - theta * v_ins * dx * dt ...
                    + 2 * kapa * T_char * theta *dx*dt / dy ...
                    / L_char / L_char);
                T_N_ins = TN(x(j));
                T_S_ins = T(i-1,j);
                if j == 1
                    aW_ins = 0;
                    aE_ins = kapa * T_char * theta * dy * dt / L2 / dx - ...
                        theta * u_ins * dy * dt / 2;
                    T_E_ins = T(i,j+1);
                    T_W_ins = TW(y(i));
                    Sp = -(dx * dy / dt - theta * v_ins * dx * dt ...
                        + 2 * kapa * T_char * theta *dx*dt / dy / L_char ...
                        / L_char + theta * u_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dy*dt ...
                        / dx / L_char / L_char);
                elseif j<M2 && j>1
                    aW_ins = kapa * T_char * theta * dy * dt / L2 / dx + ...
                        theta * u_ins * dy * dt / 2;
                    aE_ins = kapa * T_char * theta * dy * dt / L2 / dx - ...
                        theta * u_ins * dy * dt / 2;
                    T_E_ins = T(i,j+1);
                    T_W_ins = T(i,j-1);
                    Sp = -(dx * dy / dt + theta * v_ins * dx * dt ...
                        + 2 * kapa * T_char * theta *dx*dt / dy ...
                        / L_char / L_char);
                elseif j == M2
                    aW_ins = kapa * T_char * theta * dy * dt / L2 / dx + ...
                        theta * u_ins * dy * dt / 2;
                    aE_ins = 0;
                    T_E_ins = TE(y(i));
                    T_W_ins = T(i,j-1);
                    Sp = -(dx * dy / dt - theta * v_ins * dx * dt ...
                        + 2 * kapa * T_char * theta *dx*dt / dy / L_char ...
                        / L_char - theta * u_ins * dy * dt ...
                        + 2 * kapa * T_char * theta *dy*dt / dx ...
                        / L_char / L_char);
                end
            end
            Su(i,j) = aP0 * T_P_ins + aE0 * T_E_ins + aW0 * T_W_ins...
                + aN0 * T_N_ins + aS0 * T_S_ins;
            aP_ins = aN_ins + aS_ins + aW_ins + aE_ins - Sp;
            I = (j - 1) * M1 + i;
            ppp = ppp + 1;
            e(ppp) = aP_ins;
            r(ppp) = I ;
            c(ppp) = I ;
            if I<M1*M2
                ppp = ppp + 1;
                e(ppp) = -aS_ins;
                r(ppp) = I ;
                c(ppp) = I+1 ;
            end
            if I<=M1*(M2-1)
                ppp = ppp + 1;
                e(ppp) = -aE_ins;
                r(ppp) = I ;
                c(ppp) = I+M1 ;
            end
            if I>M1
                ppp = ppp + 1;
                e(ppp) = -aW_ins;
                r(ppp) = I ;
                c(ppp) = I-M1 ;
            end
            if I>1
                ppp = ppp + 1;
                e(ppp) = -aN_ins;
                r(ppp) = I ;
                c(ppp) = I-1 ;
            end
        end
    end
    matrix = sparse(r,c,e);
    re_Su = reshape_a(Su,M1,M2);
    T_downD = matrix \ re_Su;
    for j = 1 : M2
        for i = 1 : M1
            T(i,j) = T_downD((j-1)*M1+i,1);
        end
    end
    for j = 1:M2
        for i = 1:M1
            rho(i,j) = (1 - alpha * T(i,j)) * rho0(i,j);
        end
    end
    %calculate Fr,Eu,Re
    %V,L,P is characteristic
    Fr = V_char * V_char / (2 * g * L_char);
    for j = 1:M2
        for i = 1:M1
            Eu(i,j) = rho(i,j) * V_char * V_char / 2 / P_char;
            Re(i,j) = L_char * V_char * mu / rho(i,j);
        end
    end
    [u0,v0] = SIMPLE_better(Lx,Ly,M1,M2,uN,uS,vW,vE,Fr,Eu,Re,rho);
    kkk = Mt;
    if mod(kkk,10) == 0
    figure(1);
    clf;
    pcolor(x,y,T);
    hold on;
    contour(x,y,T);
    title(['2-D unsteady convection-diffusion equation @ t = ', num2str(Mt, '%.2e')])
    colorbar;
    drawnow;
    figure(2);
    quiver(x(1:3:end),y(1:3:end),u0(1:3:end,1:3:end),v0(1:3:end,1:3:end));
    startx = linspace(0,1,15);
    starty = zeros(15,1);
    streamline(x,y,u0,v0,startx,starty);
    drawnow;
    end
end