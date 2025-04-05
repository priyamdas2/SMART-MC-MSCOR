function MHat = CHat2MHat(CHat)
N = size(CHat,2);
MHat = nan(N+1, N); % Empirical transition probabilities

for i = 1:(N+1)
    MHat(i, :) = CHat(i, :) / sum(CHat(i, :));
end
end
