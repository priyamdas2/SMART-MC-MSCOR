function Total = Easom_Sum_asVec(X_vec,M)  % soln = [(1/sqrt(upto_M))*ones(M,1);zeros((M-upto_M),1)];
B = length(X_vec)/M;
X = nan(B,M);
upto_M = 2;
for b = 1:B
    X(b,:) = X_vec(((b-1)*M+1):(b*M));
end
Total = 0;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total + 1-prod(cos(sqrt(upto_M)*pi*(x(1:upto_M))))*exp(-sum((x(1:upto_M)-sqrt(1/upto_M)).^2));
end
end
    