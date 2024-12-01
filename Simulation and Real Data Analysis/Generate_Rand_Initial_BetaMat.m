% Generates (N+1)*N X (P+1) dimensonal
% matrix. Each row has L2 norm 1.
function BetaMat = Generate_Rand_Initial_BetaMat(RandSeed, N, p)
rng(RandSeed)
mu = 0;     % default           
sigma = 10; % default 
count = 0;
for i = 1:(N+1)
    for j = 1:N
        count = count + 1;
        TempArr = mu + sigma * randn(p + 1, 1);
        if(count == 1)
            BetaMat = (TempArr / sqrt(sum(TempArr.^2)))';
        else
            BetaMat = [BetaMat;(TempArr / sqrt(sum(TempArr.^2)))'];  
        end
    end
end
end