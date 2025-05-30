% Generates (N+1)*N dimensional cell each containing (P+1) dimensonal
% vector. Each vector has L2 norm 1.
function BetaCell = BetaCellRandGenerator_2(RandSeed, N, p)
rng(RandSeed)
BetaCell = cell(N+1, N);

for i = 1:(N+1)
    for j = 1:N
        TempArr = rand(p + 1, 1);
        BetaCell{i, j} = TempArr / sqrt(sum(TempArr.^2));
    end
end
end