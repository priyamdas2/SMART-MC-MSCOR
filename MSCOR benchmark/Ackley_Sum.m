function Total = Ackley_Sum(X)  % soln = (1/sqrt(M))*ones(M,1);
B = size(X,1);
M = size(X,2);
Total = 0;
for b = 1:B
    x = X(b,:)'/norm(X(b,:));
    Total = Total + (-20*exp(-0.2*sqrt(1/M*sum((x-1/sqrt(M)).^2)))-exp(1/M*sum(cos(2*pi*(x-1/sqrt(M)))))+20+exp(1));
end
end
    



