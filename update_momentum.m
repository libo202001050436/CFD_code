function [u,aP,aE,aW,aN,aS,Su] = update_momentum(Lx,Ly,...
    u,ue,uw,v,vn,vs,M,N,p,alpha,Fr,Eu,Re,type,rho,uN,uS,vE,vW,Temperature_char,beta,tempreture,rho0)

%alpha is called 亚松弛因子
%tempreture is the temperature of nodes
%beta is called the average coefficient of thermal expansion =1/V*dV/dtheta

dx = Lx / M;
dy = Ly / N;

x = dx / 2 + [0:M-1] * dx;
y = dy / 2 + [0:N-1] * dy;

aW = zeros(M,N);
aE = zeros(M,N);
aS = zeros(M,N);
aN = zeros(M,N);
aP = zeros(M,N);
Su = zeros(M,N);
re_Su = zeros(M*N,1);

Ax = dy;
Ay = dx;

ppp = 0;
r = zeros(M*N + M*(N-1)*2 + N*(M-1)*2,1);
c = zeros(M*N + M*(N-1)*2 + N*(M-1)*2,1);
e = zeros(M*N + M*(N-1)*2 + N*(M-1)*2,1);

% for i = 1:M
%     for j = 1:N
%         Dx = dy / Re(i,j) / dx;
%         Dy = dx / Re(i,j) / dy;
%     end
% end
%
% aE = zeros(M,N);
% for j = 1:N-1
%     for i = 1:M
%         aE(i,j) = Dx - dy * ue(i,j) / 2;
%     end
% end
%
% aW = zeros(M,N);
% for j = 2:N
%     for i = 1:M
%         aW(i,j) = Dx + dy * uw(i,j) / 2;
%     end
% end
%
% aN = zeros(M,N);
% for j = 1:N
%     for i = 1:M-1
%         aN(i,j) = Dy - dx * vn(i,j) / 2;
%     end
% end
%
% aS = zeros(M,N);
% for j = 1:N
%     for i = 2:M
%         aS(i,j) = Dy + dx * vs(i,j) / 2;
%     end
% end
%
% Sp = zeros(M,N);
% for j = 2:N-1
%     Sp(1,j) = -(Dy + dx * vs(1,j) / 2)*2;
%     Sp(M,j) = -(Dy - dx * vn(M,j) / 2)*2;
% end
% for i = 2:M-1
%     Sp(i,1) = -(Dx + dy * uw(i,1) / 2)*2;
%     Sp(i,N) = -(Dx + dy * ue(i,N) / 2)*2;
% end
% Sp( 1 , 1 ) = Sp(2,1)+Sp(1,2);
% Sp( 1 , N ) = Sp(2,N)+Sp(1,N-1);
% Sp( M , 1 ) = Sp(M-1,1)+Sp(M,2);
% Sp( M , N ) = Sp(M-1,N)+Sp(M,N-1);
%
% aP = aE + aW + aN + aS + dy * ue - dy * uw + dx * vn - dx * vs - Sp;

%%

for i = 1:M
    for j = 1:N
        Re_ins = Re(i,j);
        Eu_ins = Eu(i,j);

        ue_ins = ue(i,j);
        uw_ins = uw(i,j);
        un_ins = uN(x(j));
        us_ins = uS(x(j));

        vw_ins = vW(y(i));
        ve_ins = vE(y(i));
        vn_ins = vn(i,j);
        vs_ins = vs(i,j);

        u_ins = u(i,j);
        v_ins = v(i,j);

        Dx_ins = Ax / Re_ins / dx;
        Dy_ins = Ay / Re_ins / dy;

        Fe_ins = ue_ins * Ax;
        Fw_ins = uw_ins * Ax;
        Fn_ins = vn_ins * Ay;
        Fs_ins = vs_ins * Ay;

        tempreture_ins = tempreture(i,j);

        %%

            if i == 1
                aS_ins = 0;
                aN_ins = Dy_ins - Fn_ins / 2;
                if j == 1
                    aW_ins = 0;
                    aE_ins = Dx_ins - Fe_ins / 2;
                    Sp_ins = -2 * (Dy_ins + Dx_ins);
                elseif j>1 && j<N
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = Dx_ins - Fe_ins / 2;
                    Sp_ins = -2 * (Dy_ins);
                elseif j == N
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = 0;
                    Sp_ins = -2 * (Dy_ins + Dx_ins);
                end
            elseif i>1 && i<M
                aS_ins = Dy_ins + Fs_ins / 2;
                aN_ins = Dy_ins - Fn_ins / 2;
                if j == 1
                    aW_ins = 0;
                    aE_ins = Dx_ins - Fe_ins / 2;
                    Sp_ins = -2 * (Dx_ins);
                elseif j>1 && j<N
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = Dx_ins - Fe_ins / 2;
                    Sp_ins = 0;
                elseif j == N
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = 0;
                    Sp_ins = -2 * (Dx_ins);
                end
            elseif i == M
                aS_ins = Dy_ins + Fs_ins / 2;
                aN_ins = 0;
                if j == 1
                    aW_ins = 0;
                    aE_ins = Dx_ins - Fe_ins / 2;
                    Sp_ins = -2 * (Dy_ins + Dx_ins);
                elseif j>1 && j<N
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = Dx_ins - Fe_ins / 2;
                    Sp_ins = -2 * (Dy_ins);
                elseif j == N
                    aW_ins = Dx_ins + Fw_ins / 2;
                    aE_ins = 0;
                    Sp_ins = -2 * (Dy_ins + Dx_ins);
                end
            end
        
        aP_ins = aW_ins + aE_ins + aN_ins + aS_ins - Sp_ins...
            + Fe_ins - Fw_ins + Fn_ins - Fs_ins;

        I = (j - 1) * M + i;

        %record the coefficient of sparse
        ppp = ppp + 1;
        e(ppp) = aP_ins / alpha;
        r(ppp) = I ;
        c(ppp) = I ;
        if I<M*N
            ppp = ppp + 1;
            e(ppp) = -aN_ins;
            r(ppp) = I ;
            c(ppp) = I+1 ;
        end
        if I<=M*(N-1)
            ppp = ppp + 1;
            e(ppp) = -aE_ins;
            r(ppp) = I ;
            c(ppp) = I+M ;
        end
        if I>M
            ppp = ppp + 1;
            e(ppp) = -aW_ins;
            r(ppp) = I ;
            c(ppp) = I-M ;
        end
        if I>1
            ppp = ppp + 1;
            e(ppp) = -aS_ins;
            r(ppp) = I ;
            c(ppp) = I-1 ;
        end

        %for the type 1, it is the u
        if type == 1

%             u_ins = u(i,j);
            if i == 1
                if j == 1
                    Su_ins = aP_ins / alpha * (1 - alpha) * u_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dx_ins + Fw_ins) * uw_ins...
                        + (2 * Dy_ins + Fs_ins) * us_ins;
                elseif j>1 && j<N
                    Su_ins = aP_ins / alpha * (1 - alpha) * u_ins ...
                        + Ax * (p(i,j-1) - p(i,j+1)) /4 /Eu_ins...
                        + (2 * Dy_ins + Fs_ins) * us_ins;
                elseif j == N
                    Su_ins = aP_ins / alpha * (1 - alpha) * u_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dy_ins + Fs_ins) * us_ins...
                        + (2 * Dx_ins - Fe_ins) * ue_ins;
                end
            elseif i>1 && i<M
                if j == 1
                    Su_ins = aP_ins / alpha * (1 - alpha) * u_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dx_ins + Fw_ins) * uw_ins;
                elseif j>1 && j<N
                    Su_ins = aP_ins / alpha * (1 - alpha) * u_ins ...
                        + Ax * (p(i,j-1) - p(i,j+1)) /4 /Eu_ins;
                elseif j == N
                    Su_ins = aP_ins / alpha * (1 - alpha) * u_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dx_ins - Fe_ins) * ue_ins;
                end
            elseif i == M
                if j == 1
                    Su_ins = aP_ins / alpha * (1 - alpha) * u_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dy_ins - Fn_ins) * un_ins...
                        + (2 * Dx_ins + Fw_ins) * uw_ins;
                elseif j>1 && j<N
                    Su_ins = aP_ins / alpha * (1 - alpha) * u_ins ...
                        + Ax * (p(i,j-1) - p(i,j+1)) /4 /Eu_ins...
                        + (2 * Dy_ins - Fn_ins) * un_ins;
                elseif j == N
                    Su_ins = aP_ins / alpha * (1 - alpha) * u_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dy_ins - Fn_ins) * un_ins...
                        + (2 * Dx_ins - Fe_ins) * ue_ins;
                end
            end
        end

        %%

        if type == 2
%             v_ins = v(i,j);
            rho_ins = rho(i,j);
            if i == 1
                if j == 1
                    Su_ins = aP_ins / alpha * (1 - alpha) * v_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dx_ins + Fw_ins) * vw_ins...
                        + (2 * Dy_ins + Fs_ins) * vs_ins...
                        - Temperature_char * beta * tempreture_ins * (rho_ins-rho0) / 2 / Fr;
                elseif j>1 && j<N
                    Su_ins = aP_ins / alpha * (1 - alpha) * v_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dy_ins + Fs_ins) * vs_ins...
                        - Temperature_char * beta * tempreture_ins * (rho_ins-rho0) / 2 / Fr;
                elseif j == N
                    Su_ins = aP_ins / alpha * (1 - alpha) * v_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dy_ins + Fs_ins) * vs_ins...
                        + (2 * Dx_ins - Fe_ins) * ve_ins...
                        - Temperature_char * beta * tempreture_ins * (rho_ins-rho0) / 2 / Fr;
                end
            elseif i>1 && i<M
                if j == 1
                    Su_ins = aP_ins / alpha * (1 - alpha) * v_ins ...
                        + Ax * (p(i-1,j) - p(i+1,j)) /4 /Eu_ins...
                        + (2 * Dx_ins + Fw_ins) * vw_ins...
                        - Temperature_char * beta * tempreture_ins * (rho_ins-rho0) / 2 / Fr;
                elseif j>1 && j<N
                    Su_ins = aP_ins / alpha * (1 - alpha) * v_ins ...
                        + Ax * (p(i-1,j) - p(i+1,j)) /4 /Eu_ins...
                        - Temperature_char * beta * tempreture_ins * (rho_ins-rho0) / 2 / Fr;
                elseif j == N
                    Su_ins = aP_ins / alpha * (1 - alpha) * v_ins ...
                        + Ax * (p(i-1,j) - p(i+1,j)) /4 /Eu_ins...
                        + (2 * Dx_ins - Fe_ins) * ve_ins...
                        - Temperature_char * beta * tempreture_ins * (rho_ins-rho0) / 2 / Fr;
                end
            elseif i == M
                if j == 1
                    Su_ins = aP_ins / alpha * (1 - alpha) * v_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dy_ins - Fn_ins) * vn_ins...
                        + (2 * Dx_ins + Fw_ins) * vw_ins...
                        - Temperature_char * beta * tempreture_ins * (rho_ins-rho0) / 2 / Fr;
                elseif j>1 && j<N
                    Su_ins = aP_ins / alpha * (1 - alpha) * v_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dy_ins - Fn_ins) * vn_ins...
                        - Temperature_char * beta * tempreture_ins * (rho_ins-rho0) / 2 / Fr;
                elseif j == N
                    Su_ins = aP_ins / alpha * (1 - alpha) * v_ins ...
                        + Ax * (0) /2 /Eu_ins...
                        + (2 * Dy_ins - Fn_ins) * vn_ins...
                        + (2 * Dx_ins - Fe_ins) * ve_ins...
                        - Temperature_char * beta * tempreture_ins * (rho_ins-rho0) / 2 / Fr;
                end
            end
        end
    


        re_Su((j-1) * M + i,1) = Su_ins;
        Su(i,j) = Su_ins;
        aP(i,j) = aP_ins / alpha;
        aW(i,j) = aW_ins;
        aN(i,j) = aN_ins;
        aS(i,j) = aS_ins;
        aE(i,j) = aE_ins;
%         Sp(i,j) = Sp_ins;
%         Fw(i,j) = Fw_ins;
%         Fs(i,j) = Fs_ins;
%         Fn(i,j) = Fn_ins;
%         Fe(i,j) = Fe_ins;
%         Dy(i,j) = Dy_ins;
%         Dx(i,j) = Dx_ins;
    end
end
matrix = sparse(r,c,e);
% Su( 1 , 1 ) = (Su(1,2)+Su(2,1));
% Su( 1 , N ) = (Su(1,N-1)+Su(2,N));
% Su( M , 1 ) = (Su(M-1,1)+Su(M,2));
% Su( M , N ) = (Su(M-1,N)+Su(M,N-1));

% if type == 1
% re_Su = reshape_a(Su',M,N);
% end

% re_Su = reshape_a(Su',M,N);

%%

phi = matrix \ re_Su;

for j = 1:N
    for i = 1:M
        u(i,j) = phi((j-1)*M+i,1);
    end
end
end
