function [matrix] = sparse_coef_auto(aP,aW,aE,aN,aS,M,N)

re_aS = reshape_a(aS,M,N);
re_aP = reshape_a(aP,M,N);
re_aN = reshape_a(aN,M,N);
re_aW = reshape_a(aW,M,N);
re_aE = reshape_a(aE,M,N);
ppp = 0;
for i = 1:M*N
    for n = 1:5
        if n ==1 && i >= M+1
            ppp = ppp + 1;
            e(ppp) = -re_aW(i);
            r(ppp) = i ;
            c(ppp) = i-N ;
        end
        if n == 2 && i >= 2
            ppp = ppp + 1;
            e(ppp) = -re_aN(i);
            r(ppp) = i ;
            c(ppp) = i-1 ;
        end
        if n == 3
            ppp = ppp + 1;
            e(ppp) = re_aP(i);
            r(ppp) = i ;
            c(ppp) = i ;
        end
        if n == 4 && i <= M*N-1
            ppp = ppp + 1;
            e(ppp) = -re_aS(i);
            r(ppp) = i ;
            c(ppp) = 1+i ;
        end
        if n == 5 && i<=M*N-M
            ppp = ppp + 1;
            e(ppp) = -re_aE(i);
            r(ppp) = i ;
            c(ppp) = M+i ;
        end
    end
end
matrix = sparse(r,c,e);
end