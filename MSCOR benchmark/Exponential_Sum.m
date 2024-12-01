function Total = Exponential_Sum(X)   % soln = [zeros(M-1,1);1];
B = size(X,1);
M = size(X,2);
upto_M = M-1;
Total = 0;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total + 1-exp(-0.5*sum(x(1:upto_M).^2));
end
    