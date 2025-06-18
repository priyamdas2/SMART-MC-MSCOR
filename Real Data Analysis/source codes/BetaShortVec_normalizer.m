function BetaShortVec_UnitNorm = BetaShortVec_normalizer(BetaShortVecUnscaled, p)

NumNonNullLocations = length(BetaShortVecUnscaled)/(p+1);
BetaShortVec_UnitNorm = nan(NumNonNullLocations, 1);

for b = 1:NumNonNullLocations
    StartInd = (b-1)*(p+1) + 1;
    EndInd = b*(p+1);
    VecUnscaled = BetaShortVecUnscaled(StartInd:EndInd);
    BetaShortVec_UnitNorm(StartInd:EndInd) = VecUnscaled/...
        norm(VecUnscaled);
end
end