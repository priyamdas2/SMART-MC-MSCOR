function BetaMat = BetaCell2Mat(BetaCell)

N = size(BetaCell,2);
p = length(BetaCell{1,1}) - 1;
BetaMat = nan((N+1)*N,p+1);


for ii = 1:(N+1)
    for jj = 1:N
        index = (ii - 1)*N + jj;
        BetaMat(index,:) = BetaCell{ii,jj};
    end
end
end