function Total = Neg_sum_square_Sum(X)  % soln = [zeros(M-1,1);1];
B = size(X,1);
M = size(X,2);
Total = 0;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total + M-sum([1:M]*(x.^2));
end
end
    