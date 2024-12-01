function Total = Neg_log_prod_of_abs_vals_Sum(X)  % soln = (1/sqrt(M))*ones(M,1);
B = size(X,1);
M = size(X,2);
Total = 0;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total -(M/2)*log(M)-sum(log(abs(x)));
end
end
    