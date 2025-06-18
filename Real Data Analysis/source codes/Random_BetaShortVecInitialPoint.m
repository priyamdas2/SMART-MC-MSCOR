% Based on provided Random Generator Seed Number, and the indicator matrix
% of non-NULL locations, it creates BetaShortVec of dimension 
% = (number of non-NULL locations)*(p+1)

function BetaShortVecInitial = Random_BetaShortVecInitialPoint(RandSeed, NonNullPositions, p)
rng(RandSeed)
N = size(NonNullPositions, 2);
BetaCell = cell(N+1, N);
mu = 0;
sigma = 1;

for i = 1:(N+1)
    for j = 1:N
        TempArr = mu + sigma * randn(p + 1, 1);
        BetaCell{i, j} = TempArr / sqrt(sum(TempArr.^2));
    end
end

BetaVec = Beta_cell2vec(BetaCell);
BetaShortVecInitial = BetaVec2ShortVec(BetaVec,NonNullPositions);
end