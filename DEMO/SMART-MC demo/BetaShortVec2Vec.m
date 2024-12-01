% Takes a (num_non_null_positions) * (p+1)  dinemsinal beta vector, and 
% reconstucts the full vector by incorporating 0 at NULL locations. 
% Length of final beta is (N+1)*N*(p+1)
function BetaVecFinal = BetaShortVec2Vec(BetaShortVec,NonNullPositions)

N = size(NonNullPositions, 2);
NumNonZero = sum(sum(NonNullPositions));
p = length(BetaShortVec)/NumNonZero - 1;

count = 0;
BetaVec = nan;
for i = 1:(N+1)
    for j = 1:N
        if(NonNullPositions(i,j) == 1)
            StartInd = count + 1;
            EndInd = count + p + 1;
            BetaVec = [BetaVec; BetaShortVec(StartInd:EndInd)];
            count = count + p + 1;
        else
            BetaVec = [BetaVec; zeros(p+1,1)];
        end
    end
end

BetaVecFinal = BetaVec;
BetaVecFinal(1) = [];

end


