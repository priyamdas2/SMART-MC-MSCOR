function Total = Rastrigin_Sum_asVec(X_vec,M)  % soln = -(1/sqrt(M))*ones(M,1);
B = length(X_vec)/M;
X = nan(B,M);
for b = 1:B
    X(b,:) = X_vec(((b-1)*M+1):(b*M));
end
Total = 0;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total + (10*M+sum((x+1/sqrt(M)).^2-10*cos(2*pi*(x+1/sqrt(M)))));
end
end
    



