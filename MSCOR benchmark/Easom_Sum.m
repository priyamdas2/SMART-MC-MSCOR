function Total = Easom_Sum(X)  % soln = [(1/sqrt(upto_M))*ones(M,1);zeros((M-upto_M),1)];
B = size(X,1);
M = size(X,2);
upto_M = 2;
Total = 0;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total + 1-prod(cos(sqrt(upto_M)*pi*(x(1:upto_M))))*exp(-sum((x(1:upto_M)-sqrt(1/upto_M)).^2));
end
end
    