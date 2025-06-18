function [CoeffsEst, CoeffsTrue,DiffMatTrans, RMSE, MAD]=...
    Simulation_TopHowManyBetaComparison(TopHowManyInitials,...
    TopHowManyTrans, SeqCell, BetaCellEst,BetaCellTrue,InterceptIndic, N) 
p = max(size(BetaCellEst{1,1})) - 1;
CHat = TransitionCountMat(SeqCell, N);

%%% Initial State Probabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHatInitial = CHat(1,:)';
NumTotalInitials = sum(CHatInitial);
[ValuesIn,ranksIn] = sort(CHatInitial,'descend');
TopIndicesIn = ranksIn(1:TopHowManyInitials);


CoeffsInitialsEst = nan(TopHowManyInitials,4+p+1);
CoeffsInitialsTrue = nan(TopHowManyInitials,4+p+1);

for i = 1:TopHowManyInitials
    CoeffsInitialsEst(i,1) = TopIndicesIn(i);
    CoeffsInitialsEst(i,3) = CHatInitial(TopIndicesIn(i));
    CoeffsInitialsEst(i,4) = 100*CoeffsInitialsEst(i,3)/NumTotalInitials; % percentage
    CoeffsInitialsEst(i,5:end) = BetaCellEst{1, TopIndicesIn(i)}';
    
    CoeffsInitialsTrue(i,1:4) = CoeffsInitialsEst(i,1:4);
    CoeffsInitialsTrue(i,5:end) = BetaCellTrue{1, TopIndicesIn(i)}';
end
    
%%% Transition Probabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHatTrans = CHat(2:end,:);
NumTotalTransitions = sum(sum(CHatTrans));
[Values,ranks] = sort(CHatTrans(:),'descend');
TopIndices = ranks(1:TopHowManyTrans);
TopColInd = floor((TopIndices-0.00001)/N) + 1;
TopRowInd = TopIndices - N * floor((TopIndices-0.00001)/N);
% TreatmentNamesALL = ["1","2","3","4","5","6","7","8","9","10","11","12",...
%    "13","14","15","16","17","18","19","20"];
% TreatmentNamesShort = TreatmentNamesALL(1:N);


CoeffsTransEst = nan(TopHowManyTrans,4+p+1);
CoeffsTransTrue = nan(TopHowManyTrans,4+p+1);
for i = 1:TopHowManyTrans
    CoeffsTransEst(i,1:2) = [TopRowInd(i), TopColInd(i)];
    CoeffsTransEst(i,3) = CHatTrans(TopRowInd(i), TopColInd(i));
    CoeffsTransEst(i,4) = 100*CoeffsTransEst(i,3)/NumTotalTransitions; % percentage
    CoeffsTransEst(i,5:end) = BetaCellEst{TopRowInd(i)+1, TopColInd(i)}';
    
    CoeffsTransTrue(i,1:4) = CoeffsTransEst(i,1:4);
    CoeffsTransTrue(i,5:end) = BetaCellTrue{TopRowInd(i)+1, TopColInd(i)}';
end
    

CoeffsEstTemp = [CoeffsInitialsEst;CoeffsTransEst];
CoeffsTrueTemp = [CoeffsInitialsTrue;CoeffsTransTrue];

if(InterceptIndic == 0)
    CoeffsEstTemp2 = CoeffsEstTemp;
    CoeffsTrueTemp2 = CoeffsTrueTemp;
    CoeffsEstTemp2(:,5) = [];
    CoeffsTrueTemp2(:,5) = [];
    CoeffsEst = CoeffsEstTemp2;
    CoeffsTrue = CoeffsTrueTemp2;
    NumRows = TopHowManyInitials+TopHowManyTrans; 
    for ii = 1:NumRows
        CoeffsEst(ii,5:end) = CoeffsEstTemp2(ii,5:end)/norm(CoeffsEstTemp2(ii,5:end));
        CoeffsTrue(ii,5:end) = CoeffsTrueTemp2(ii,5:end)/norm(CoeffsTrueTemp2(ii,5:end));
    end
else
    CoeffsEst = CoeffsEstTemp;
    CoeffsTrue = CoeffsTrueTemp;
end
        

DiffMatTrans = abs(CoeffsEst(2:end,5:end) - CoeffsTrue(2:end,5:end)); % Excluding Initial state vector diffs
RMSE = sqrt(sum(sum(DiffMatTrans.^2))/(size(DiffMatTrans,1)*size(DiffMatTrans,2)));
MAD = sum(sum(DiffMatTrans))/(size(DiffMatTrans,1)*size(DiffMatTrans,2));

    
end
