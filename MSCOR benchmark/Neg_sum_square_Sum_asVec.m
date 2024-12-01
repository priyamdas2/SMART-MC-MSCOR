function Total = Neg_sum_square_Sum_asVec(X_vec,M)  % soln = [zeros(M-1,1);1];
B = length(X_vec)/M;
X = nan(B,M);
for b = 1:B
    X(b,:) = X_vec(((b-1)*M+1):(b*M));
end
Total = 0;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total + M-sum([1:M]*(x.^2));
end
end
    