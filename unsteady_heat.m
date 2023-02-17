function [theta,u,v] = unsteady_heat(Lx,Ly,time,dt,T0,TW,TE,TN,TS,dx,dy,rho_background,...
    alpha,V_char,L_char,T_char,mu,g,P_char,Temperature_char,heat_cap,k_inter,uN,uS,vW,vE,WP)

%k is the coefficient of thermal conductive
%theta is the temperature of nodes
%alpha is the average coefficient of thermal expansion =1/V*dV/dtheta
%T_char relys on time not temperature
ppp = 0;

M1 = Lx / dx;
x  = dx / 2 + [0:M1-1] * dx;

M2 = Ly / dy;
y  = dy / 2 + [0:M2-1] * dy;

Mt = time / dt;

u = zeros(M1,M2);
v = zeros(M1,M2);
theta = zeros(M1,M2);
k = zeros(M1,M2);
for i = 1:M1
    for j = 1:M2
        theta(i,j) = T0;
%         if i > 23 && i < 27 && j > 23 && j < 27
%             theta(i,j) = 1;
%         end
    end
end
for i = 1:M1
    for j = 1:M2
        %                 if j == 1 || j == M2
        %                     k(i,j) = 1e-9;
        %                 elseif 1
        k(i,j) = k_inter;
        %                 end
    end
end


Pe = zeros(M1,M2);
St = zeros(M1,M2);
rho= zeros(M1,M2);
Eu = zeros(M1,M2);
Re = zeros(M1,M2);

aP = zeros(M1,M2);
aW = zeros(M1,M2);
aE = zeros(M1,M2);
aN = zeros(M1,M2);
aS = zeros(M1,M2);
Su = zeros(M1,M2);

nu = zeros(M1,M2);



for j = 1:M2
    for i = 1:M1
        rho(i,j) = rho_background;
    end
end

for t = 0:Mt

    %give the initial value of v, u, T
    if t > 0
        u = u0;
        v = v0;
    end
    for i = 1:M1
        for j = 1:M2
            %                 theta(i,1) = TW(y(i));
            %                 theta(i,M2) = TE(y(i));
            %                 theta(1,j) = TS(y(i));
            %                 theta(M1,j) = TN(y(i));
            if j == 1
                v(i,j) = vW(y(i));
                if i == 1
                    v(i,j) = 0;
                elseif i == M1
                    v(i,j) = 0;
                end
            elseif j == M2
                v(i,j) = vE(y(i));
                if i == 1
                    v(i,j) = 0;
                elseif i == M1
                    v(i,j) = 0;
                end
            end
            if i == 1
                u(i,j) = uS(x(j));
                if j == 1
                    u(i,j) = 0;
                elseif j == M2
                    u(i,j) = 0;
                end
            elseif i == M1
                u(i,j) = uN(x(j));
                if j == 1
                    u(i,j) = 0;
                elseif j == M2
                    u(i,j) = 0;
                end
            end
        end
    end

    %%give the value to the contemporary field:u0, v0, theta0
    u0 = u;
    v0 = v;
    theta0 = theta;
    rho0 = rho;


    for j = 1:M2
        for i = 1:M1
            rho(i,j) = (1 - alpha * theta0(i,j)) * rho0(i,j);
        end
    end
    for j = 1:M2
        for i = 1:M1
            nu(i,j) = mu / rho(i,j);
        end
    end
    for j = 1:M2
        for i = 1:M1
            Eu(i,j) = rho(i,j) * V_char * V_char / 2 / P_char;
            Re(i,j) = L_char * L_char * nu(i,j) / rho(i,j) / T_char;
        end
    end
    Fr = V_char * V_char / 2 / g / L_char;



    %%calculate the new value of v, u
    [u,v,ue,uw,vn,vs] = SIMPLE_better(Lx,Ly,M1,M2,uN,uS,vW,vE,Fr,Eu,Re,rho,Temperature_char,alpha,theta0,rho_background);


    %%calculate the value of coefficient
    for i = 1:M1
        for j = 1:M2
            Pe(i,j) = L_char * V_char * rho(i,j) * heat_cap / k(i,j);
            St(i,j) = L_char / V_char / T_char;
        end
    end


    if t > 0
        aE0 = aE;
        aW0 = aW;
        aN0 = aN;
        aS0 = aS;
        aP0 = aP;
    end

    %%calculate the temperature coefficients of metrix
    for i = 1:M1
        for j = 1:M2

            %% new value

            ue_ins = ue(i,j);
            uw_ins = uw(i,j);
            vn_ins = vn(i,j);
            vs_ins = vs(i,j);

            theta0_P_ins = theta0(i,j);

            Pe_ins = Pe(i,j);
            St_ins = St(i,j);


            Dx_ins = dy / Pe_ins / dx;
            Dy_ins = dx / Pe_ins / dy;

            Fe_ins = ue_ins * dy;
            Fw_ins = uw_ins * dy;
            Fn_ins = vn_ins * dx;
            Fs_ins = vs_ins * dx;

            %% calculate temperature

            if i == 1
                aS_ins = 0;
                aN_ins = Dy_ins - Fn_ins / 2;

                theta0_N_ins = theta0(i+1,j);
                theta0_S_ins = TS(x(j));
                if j == 1
                    aW_ins = Fw_ins / 2;
                    aE_ins = Dx_ins - Fe_ins / 2;

                    theta0_E_ins = theta0(i,j+1);
                    theta0_W_ins = theta0(i,j) * 2 - theta0(i,j+1);

                    Sp_ins = -2 * (Dy_ins + Dx_ins);

                    variable = (2 * Dx_ins) * theta0_W_ins...
                        + (2 * Dy_ins) * theta0_S_ins;
                elseif j>1 && j<M2
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = Dx_ins - Fe_ins / 2;

                    theta0_E_ins = theta0(i,j+1);
                    theta0_W_ins = theta0(i,j-1);

                    Sp_ins = -2 * (Dy_ins);

                    variable = (2 * Dy_ins) * theta0_S_ins;
                elseif j == M2
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = - Fe_ins / 2;

                    theta0_E_ins = theta0(i,j) * 2 - theta0(i,j-1);
                    theta0_W_ins = theta0(i,j-1);

                    Sp_ins = -2 * (Dy_ins + Dx_ins);

                    variable = (2 * Dx_ins) * theta0_E_ins...
                        + (2 * Dy_ins) * theta0_S_ins;
                end
            elseif i>1 && i<M1
                aS_ins = Dy_ins + Fs_ins / 2;
                aN_ins = Dy_ins - Fn_ins / 2;

                theta0_N_ins = theta0(i+1,j);
                theta0_S_ins = theta0(i-1,j);
                if j == 1
                    aW_ins = Fw_ins / 2;
                    aE_ins = Dx_ins - Fe_ins / 2;

                    theta0_E_ins = theta0(i,j+1);
                    theta0_W_ins = theta0(i,j) * 2 - theta0(i,j+1);

                    Sp_ins = -2 * (0);

                    variable = (2 * 0) * theta0_W_ins;
                elseif j>1 && j<M2
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = Dx_ins - Fe_ins / 2;

                    theta0_E_ins = theta0(i,j+1);
                    theta0_W_ins = theta0(i,j-1);

                    Sp_ins = 0;

                    variable = 0;
                elseif j == M2
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = 0 - Fe_ins / 2;

                    theta0_E_ins = theta0(i,j) * 2 - theta0(i,j-1);
                    theta0_W_ins = theta0(i,j-1);

                    Sp_ins = 0;

                    variable = (2 * 0) * theta0_E_ins;
                end
            elseif i == M1
                aS_ins = Dy_ins + Fs_ins / 2;
                aN_ins = 0;

                theta0_N_ins = theta0(i,j) * 2 - theta0(i-1,j);
                theta0_S_ins = theta0(i-1,j);
                if j == 1
                    aW_ins = Fw_ins / 2;
                    aE_ins = Dx_ins - Fe_ins / 2;

                    theta0_E_ins = theta0(i,j+1);
                    theta0_W_ins = theta0(i,j) * 2 - theta0(i,j+1);

                    Sp_ins = -2 * (Dy_ins);

                    variable = (2 * 0) * theta0_W_ins...
                        + (2 * Dy_ins - Fn_ins) * theta0_N_ins;
                elseif j>1 && j<M2
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = Dx_ins - Fe_ins / 2;

                    theta0_E_ins = theta0(i,j+1);
                    theta0_W_ins = theta0(i,j-1);

                    Sp_ins = -2 * (Dy_ins);

                    variable = (2 * Dy_ins) * theta0_N_ins;
                elseif j == M2
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = - Fe_ins / 2;

                    theta0_E_ins = theta0(i,j) * 2 - theta0(i,j-1);
                    theta0_W_ins = theta0(i,j-1);

                    Sp_ins = -2 * (Dy_ins + 0);

                    variable = (2 * 0) * theta0_E_ins...
                        + (2 * Dy_ins) * theta0_N_ins;
                end
            end

            aP_ins = aW_ins + aE_ins + aN_ins + aS_ins - Sp_ins...
                + Fe_ins - Fw_ins + Fn_ins - Fs_ins + St_ins * dt * dx * dy;
            aP(i,j) = aP_ins;
            aN(i,j) = aN_ins;
            aS(i,j) = aS_ins;
            aW(i,j) = aW_ins;
            aE(i,j) = aE_ins;

            if t == 0
                ue_ins = u0(i,j);
                uw_ins = u0(i,j);
                vn_ins = v0(i,j);
                vs_ins = v0(i,j);

                theta0_P_ins = theta0(i,j);

                Pe_ins = Pe(i,j);
                St_ins = St(i,j);


                Dx_ins = dy / Pe_ins / dx;
                Dy_ins = dx / Pe_ins / dy;

                Fe_ins = ue_ins * dy;
                Fw_ins = uw_ins * dy;
                Fn_ins = vn_ins * dx;
                Fs_ins = vs_ins * dx;

                %% calculate temperature

                if i == 1
                    aS0_ins = 0;
                    aN0_ins = Dy_ins - Fn_ins / 2;
                    if j == 1
                        aW0_ins = 0;
                        aE0_ins = 0 - Fe_ins / 2;
                        Sp0_ins = -2 * (Dy_ins + 0);
                    elseif j>1 && j<M2
                        aW0_ins = Dx_ins + Fw_ins / 2;
                        aE0_ins = Dx_ins - Fe_ins / 2;
                        Sp0_ins = -2 * (Dy_ins);
                    elseif j == M2
                        aW0_ins = 0 + Fw_ins / 2;
                        aE0_ins = 0;
                        Sp0_ins = -2 * (Dy_ins + 0);
                    end
                elseif i>1 && i<M1
                    aS0_ins = Dy_ins + Fs_ins / 2;
                    aN0_ins = Dy_ins - Fn_ins / 2;
                    if j == 1
                        aW0_ins = 0;
                        aE0_ins = 0 - Fe_ins / 2;
                        Sp0_ins = -2 * (0);
                    elseif j>1 && j<M2
                        aW0_ins = Dx_ins + Fw_ins / 2;
                        aE0_ins = Dx_ins - Fe_ins / 2;
                        Sp0_ins = 0;
                    elseif j == M2
                        aW0_ins = 0 + Fw_ins / 2;
                        aE0_ins = 0;
                        Sp0_ins = -2 * (0);
                    end
                elseif i == M1
                    aS_ins = Dy_ins + Fs_ins / 2;
                    aN_ins = 0;
                    if j == 1
                        aW0_ins = 0;
                        aE0_ins = 0 - Fe_ins / 2;
                        Sp0_ins = -2 * (Dy_ins + 0);
                    elseif j>1 && j<M2
                        aW0_ins = Dx_ins + Fw_ins / 2;
                        aE0_ins = Dx_ins - Fe_ins / 2;
                        Sp0_ins = -2 * (Dy_ins);
                    elseif j == M2
                        aW0_ins = 0 + Fw_ins / 2;
                        aE0_ins = 0;
                        Sp0_ins = -2 * (Dy_ins + 0);
                    end
                end
            aP0_ins = aW0_ins + aE0_ins + aN0_ins + aS0_ins - Sp0_ins...
                + Fe_ins - Fw_ins + Fn_ins - Fs_ins + St_ins * dt * dx * dy;
                timegoing = WP;
            elseif t
                aE0_ins = aE0(i,j);
                aW0_ins = aW0(i,j);
                aN0_ins = aN0(i,j);
                aS0_ins = aS0(i,j);
                aP0_ins = aP0(i,j);
                timegoing = 1;
            end


            Su_ins = (1 - WP) / WP * (aE0_ins * theta0_E_ins +...
                aW0_ins * theta0_W_ins + aN0_ins * theta0_N_ins +...
                aS0_ins * theta0_S_ins - aP0_ins * theta0_P_ins) +...
                variable / WP * timegoing;

            %             if j == 1 || j == M2 || i == M1
            %                 Su_ins = 0;
            %             end
            varys(i,j) = variable;
            theta0_W(i,j) = theta0_W_ins;
            theta0_E(i,j) = theta0_E_ins;
            theta0_N(i,j) = theta0_N_ins;
            theta0_S(i,j) = theta0_S_ins;
            theta0_P(i,j) = theta0_P_ins;





            Su(i,j) = Su_ins;


            I = (j - 1) * M1 + i;

            %record the coefficient of sparse
            ppp = ppp + 1;
            e(ppp) = aP_ins;
            r(ppp) = I ;
            c(ppp) = I ;
            if I<M1*M2
                ppp = ppp + 1;
                e(ppp) = -aN_ins;
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
                e(ppp) = -aS_ins;
                r(ppp) = I ;
                c(ppp) = I-1 ;
            end

            re_Su((j-1) * M1 + i,1) = Su_ins;

        end
    end

    matrix = sparse(r,c,e);

    phi = matrix \ re_Su;

    for i = 1:M1
        for j = 1:M2
            theta(i,j) = phi((j-1)*M1+i,1);
        end
    end

    kkk = t;
    if mod(kkk,1) == 0
%                         figure(1);
%                         subplot(3,3,1);
%                         mesh(y,x,theta0_E);
%                         xlabel('x(m)');
%                         ylabel('y(m)');
%                         zlabel('theta0 E');
%                         box on;
%                         figure(1);
%                         subplot(3,3,2);
%                         mesh(y,x,theta0_W);
%                         xlabel('x(m)');
%                         ylabel('y(m)');
%                         zlabel('theta0 W');
%                         box on;
%                         figure(1);
%                         subplot(3,3,3);
%                         mesh(y,x,theta0_N);
%                         xlabel('x(m)');
%                         ylabel('y(m)');
%                         zlabel('theta0 N');
%                         box on;
%                         figure(1);
%                         subplot(3,3,4);
%                         mesh(y,x,theta0_S);
%                         xlabel('x(m)');
%                         ylabel('y(m)');
%                         zlabel('theta0 S');
%                         box on;
%                         figure(1);
%                         subplot(3,3,5);
%                         mesh(y,x,theta0_P);
%                         xlabel('x(m)');
%                         ylabel('y(m)');
%                         zlabel('theta0 P');
%                         box on;
%                         figure(1);
%                         subplot(3,3,6);
%                         mesh(y,x,varys);
%                         xlabel('x(m)');
%                         ylabel('y(m)');
%                         zlabel('varys');
%                         box on;
%                         figure(1);
%                         subplot(3,3,7);
%                         mesh(y,x,aP);
%                         xlabel('x(m)');
%                         ylabel('y(m)');
%                         zlabel('aP');
%                         box on;
%                         figure(1);
%                         subplot(3,3,8);
%                         mesh(y,x,theta);
%                         xlabel('x(m)');
%                         ylabel('y(m)');
%                         zlabel('theta');
%                         box on;
%                         figure(1);
%                         subplot(3,3,9);
%                         mesh(y,x,Su);
%                         xlabel('x(m)');
%                         ylabel('y(m)');
%                         zlabel('Su');
%                         box on;
        clf;
        %         figure(1);
        %         subplot(5,2,10);
        % pcolor(x,y,theta);
        
        
        % hold on;
        contour(x,y,theta, 'LineColor','k','LevelStep',0.05,'Fill','on');
        %         title(['2-D unsteady convection-diffusion equation step t = ', num2str(t, '%d')]);
        shading interp;
        title(['step t = ', num2str(t, '%d')]);
        colorbar;
        hold on;
        quiver(x(1:3:end),y(1:3:end),u(1:3:end,1:3:end),v(1:3:end,1:3:end));
        hold on;
        drawnow;
        %     startx = linspace(0,3,15);
        %     starty = startx;
        %     streamline(x,y,u,v,startx,starty);
        %     hold on
        str1 = ['NSH的第',num2str(t),'个图'];
        print(gcf,[str1,'.png'],'-dpng');
    end
end
end
