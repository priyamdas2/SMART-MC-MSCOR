% Converts Beta cell to an (N+1)*N*(p+1) diemnsional array,
% traversing ROW-wise across cell blocks
function BetaVecFinal = Beta_cell2vec(BetaCell)
N = size(BetaCell, 2);
p = size(BetaCell{1,1}, 1) - 1;
BetaVec = nan;
for i = 1:(N+1)
    for j = 1:N
        BetaVec = [BetaVec; BetaCell{i,j}];
    end
end
BetaVecFinal = BetaVec;
BetaVecFinal(1) = [];
end
