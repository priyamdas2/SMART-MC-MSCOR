function [CoeffSummary, SeSummary] = TopHowManyTransFormula(TopHowMany, SeqCell, BetaCellFinal, BetaALLBootStdErr, N)

CHat = TransitionCountMat(SeqCell, N);
CHatTrans = CHat(2:end,:);


[Values,ranks] = sort(CHatTrans(:),'descend');
TopIndices = ranks(1:TopHowMany);
TopColInd = floor((TopIndices-0.00001)/N) + 1;
TopRowInd = TopIndices - N * floor((TopIndices-0.00001)/N);
TreatmentNamesShort = ["BcD", "GA", "IB","DF","Nat","Fin","Ter","Cyc",...
    "Mit","Ale"];


for i = 1:TopHowMany
    [TopRowInd(i), TopColInd(i)];
    CHatTrans(TopRowInd(i), TopColInd(i));
    fprintf('Rank %d: Transition from %s to %s (Transition count = %d) \n',...
        i, TreatmentNamesShort(TopRowInd(i)),...
        TreatmentNamesShort(TopColInd(i)),...
        CHatTrans(TopRowInd(i), TopColInd(i)));
    fprintf('Covariate dependency: %.2f* Int + %.2f * ageT + %.2f * dis.dur.T + %.2f * genderF + %.2f * white + %.2f * black \n',...
        BetaCellFinal{TopRowInd(i)+1, TopColInd(i)}); % +1 since we use only transintion mat
    
    %%% Extract Bootstrap Std. Errors
    ind_row = TopRowInd(i);
    ind_col = TopColInd(i);
    len = size(BetaALLBootStdErr,3);
    StdErrs = nan(len,1);
    for ii = 1:len
        StdErrs(ii) = BetaALLBootStdErr(ind_row + 1, ind_col, ii); % +1 since we use only transintion mat
    end
    fprintf('Bootstrap Standard errors are: %.3f (Int), %.3f (ageT), %.3f (dis.dur.T), %.3f (genderF), %.3f (white), %.3f (black). \n',...
        StdErrs); % +1 since we use only transintion mat
    fprintf('\n');
    
    if(i == 1)
        CoeffSummary = [TreatmentNamesShort(TopRowInd(i)),...
            TreatmentNamesShort(TopColInd(i)),...
            CHatTrans(TopRowInd(i), TopColInd(i)),...
            round(BetaCellFinal{TopRowInd(i)+1, TopColInd(i)}',2)];
        SeSummary = [TreatmentNamesShort(TopRowInd(i)),...
            TreatmentNamesShort(TopColInd(i)),...
            CHatTrans(TopRowInd(i), TopColInd(i)),round(StdErrs',3)];
        
    else
        CoeffSummary = [CoeffSummary;TreatmentNamesShort(TopRowInd(i)),...
            TreatmentNamesShort(TopColInd(i)),...
            CHatTrans(TopRowInd(i), TopColInd(i)),...
            round(BetaCellFinal{TopRowInd(i)+1, TopColInd(i)}',2)];
        
        SeSummary = [SeSummary; TreatmentNamesShort(TopRowInd(i)),...
            TreatmentNamesShort(TopColInd(i)),...
            CHatTrans(TopRowInd(i), TopColInd(i)),round(StdErrs',3)];
    end
end
end
