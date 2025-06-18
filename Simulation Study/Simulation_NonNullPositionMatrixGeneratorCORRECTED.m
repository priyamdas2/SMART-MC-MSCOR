function NonNullPositionsTrue = Simulation_NonNullPositionMatrixGeneratorCORRECTED(RandSeed,TrueNullPercentage,N)

rng(RandSeed)
VecDiagIndices = nan(N,1);
VecOneNondiagIndices = nan(N,1);
for j =1:N
    VecDiagIndices(j) = (N+1)*(j-1)+1+j;                                   % Making sure diag transition probs are non-zero
    RestColIndexes = setdiff(1:N,j);
    RandNondiagIndex = randsample(RestColIndexes,1);
    VecOneNondiagIndices(j) = (N+1)*(RandNondiagIndex-1)+1+j;              % Making sure one non-diag transition probs are non-zero
end
NonNullInitials = randsample(N,2,false);
VecTwoNonNullInitials = (NonNullInitials-1)*(N+1)+1;                       % Making sure two initial probs are non-zero
VecDontDeleteIndices = [VecTwoNonNullInitials;VecDiagIndices;VecOneNondiagIndices];


VecIndices = 1:((N+1)*N);
VecCanBeDeletedIndices = setdiff(VecIndices,VecDontDeleteIndices);
NumIndices2Delete = round((TrueNullPercentage/100)*(N+1)*N);
Indices2Delete = randsample(VecCanBeDeletedIndices, NumIndices2Delete, false);
NonNullPositionsTrue = ones(N+1,N);
NonNullPositionsTrue(Indices2Delete) = 0;
end