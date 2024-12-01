function Total = Griewank_Sum(X)   % soln = (1/sqrt(M))*ones(M,1);
B = size(X,1);
M = size(X,2);
Total = 0;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total +(1/4000)*sum((x-1/sqrt(M)).^2) - prod(cos(transpose((x-1/sqrt(M)))./sqrt(1:M)))+1;
end
end
    



