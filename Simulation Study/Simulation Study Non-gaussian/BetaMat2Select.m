function BetaMatSelected = BetaMat2Select(BetaMat,N, NonNullPositions)
p = size(BetaMat,2) - 1;
BetaMatSelected = nan(sum(sum(NonNullPositions)),p+1);

tempIdx = 0;
for rIdx = 1:(N+1)
    for cIdx = 1:N
        serialNum = (rIdx-1)*N + cIdx;
        if(NonNullPositions(rIdx,cIdx) == 1)
            tempIdx = tempIdx+1;
            if(tempIdx == 1)
                BetaMatSelected = BetaMat(serialNum,:);
            else
                BetaMatSelected = [BetaMatSelected; BetaMat(serialNum,:)];
            end
        end
    end
end
end