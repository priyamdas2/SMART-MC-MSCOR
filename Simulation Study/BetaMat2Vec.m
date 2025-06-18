function BetaVec = BetaMat2Vec(BetaMat)
NumRows = size(BetaMat,1);

for i = 1:NumRows
    if(i == 1)
        BetaVec = BetaMat(i,:)';
    else
        BetaVec = [BetaVec; BetaMat(i,:)'];
    end
end
end