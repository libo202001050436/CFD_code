function [re_Su] = reshape_a(Su,M,N)

re_Su = zeros(M*N,1);
for j = 1 : N
    for i = 1 : M
        re_Su((j-1)*M+i,1) = Su(i,j);
    end
end
return