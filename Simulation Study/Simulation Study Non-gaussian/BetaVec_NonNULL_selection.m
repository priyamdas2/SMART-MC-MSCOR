% Takes Beta vector of length (N+1)*N*(p+1), DELETES the values
% corresponding to NULL location, then returns the rest
function BetaVec = BetaVec_NonNULL_selection(BetaVecWhole, NonNullPositions)

N = size(NonNullPositions, 2);
p = length(BetaVecWhole)/((N+1)*N) - 1;

BetaCellTemp = Beta_vec2cell(BetaVecWhole, N);

for i = 1:(N+1)
    for j = 1:N
        if(NonNullPositions(i, j) == 0)
            BetaCellTemp{i,j} = zeros(p+1, 1);
        end
    end
end
BetaVec = Beta_cell2vec(BetaCellTemp);
end
