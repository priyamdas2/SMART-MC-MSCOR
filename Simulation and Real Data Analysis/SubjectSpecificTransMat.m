function SubjectTransMats = SubjectSpecificTransMat(X,BetaCell,NonNullPositionsTrue)
NumSubjects = size(X,1);
N =  size(BetaCell,2);
SubjectTransMats = nan(N+1,N,NumSubjects);
TempMat = nan(N+1,N,NumSubjects);
for k = 1:NumSubjects
    for i = 1:(N+1)
        for j = 1:N
            TempMat(i,j,k) = exp(X(k,:)*BetaCell{i,j});
            if(NonNullPositionsTrue(i,j) == 0)
                TempMat(i,j,k) = 0;
            end
        end
        SubjectTransMats(i,:,k) = TempMat(i,:,k)./sum(TempMat(i,:,k));
    end
end
end
