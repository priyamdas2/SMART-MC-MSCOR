% Takes a (N+1)*N*(p+1) dinemsinal vector, and retains the vector values
% corresponding to the NON-NULL locations. Length of short_vec is
% = num_non_null_positions * (p+1)
function BetaShortVecFinal = BetaVec2ShortVec(BetaVec,NonNullPositions)
N = size(NonNullPositions, 2);
p = length(BetaVec)/((N+1)*N) - 1;
BetaShortVec = nan;
for i = 1:(N+1)
    for j = 1:N
        if(NonNullPositions(i,j) == 1)
            StartInd = (i-1)*N*(p+1)+(j-1)*(p+1) + 1;
            EndInd = (i-1)*N*(p+1)+(j-1)*(p+1) + p + 1;
            BetaShortVec = [BetaShortVec; BetaVec(StartInd:EndInd)];
        end
    end
end

BetaShortVecFinal = BetaShortVec;
BetaShortVecFinal(1) = [];

if(sum(sum(NonNullPositions)) ~= (length(BetaShortVecFinal) / (p+1)))
    print("Error");
end
end
        