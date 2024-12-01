function Total = Exponential_Sum_asVec(X_vec,M)   % soln = [zeros(M-1,1);1];
B = length(X_vec)/M;
X = nan(B,M);
for b = 1:B
    X(b,:) = X_vec(((b-1)*M+1):(b*M));
end
Total = 0;
upto_M = M-1;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total + 1-exp(-0.5*sum(x(1:upto_M).^2));
end
    