function TrtSeqMat = Simulation_generate_TrtSeqMat(RandSeed,SubjectTransMats,TrtSeqLengthEach)

rng(RandSeed)
NumSubjects = size(SubjectTransMats,3);
TrtSeqMat = nan(NumSubjects*TrtSeqLengthEach,2);
for k = 1:NumSubjects
    IntialTrt = find(mnrnd(1, SubjectTransMats(1,:,k), 1) == 1);
    LastTrt = IntialTrt;
    CumTrt = IntialTrt;
    for t = 2:TrtSeqLengthEach
        LastTrt = CumTrt(end);
        ThisTrt = find(mnrnd(1, SubjectTransMats(LastTrt+1,:,k), 1) == 1);
        CumTrt = [CumTrt;ThisTrt];
    end
    startIdx = (k-1)*TrtSeqLengthEach + 1;
    endIdx = k*TrtSeqLengthEach;
    TrtSeqMat(startIdx:endIdx,1) = k;
    TrtSeqMat(startIdx:endIdx,2) = CumTrt;
end
end