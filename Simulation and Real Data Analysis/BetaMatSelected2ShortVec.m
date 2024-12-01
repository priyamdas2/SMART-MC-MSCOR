function BetaShortVec = BetaMatSelected2ShortVec(BetaMatSelected)
NumNonNullPositions = size(BetaMatSelected,1);

for ii = 1:NumNonNullPositions
    if(ii == 1)
        BetaShortVec = BetaMatSelected(ii,:)';
    else
        BetaShortVec =[BetaShortVec; BetaMatSelected(ii,:)'];
    end
end