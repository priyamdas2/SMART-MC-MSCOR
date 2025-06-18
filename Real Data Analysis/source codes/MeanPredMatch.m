function MeanPredMatchOverall = MeanPredMatch(M_opt, SeqCell, RandSeed, NumExp)
rng(RandSeed)
K = size(SeqCell,2); % Number of patients

NumTransitions = zeros(K,1);
MatchPropSum = zeros(K,1);

for k = 1:K
    TrueStates = SeqCell{k};
    t_k = length(TrueStates);
    % ProportionMatch
    for t = 1:t_k
        if(t == 1)
            Obs = mnrnd(NumExp,M_opt(1,:,k));
            MatchProp = Obs(TrueStates(t))/NumExp;
        else
            Obs = mnrnd(NumExp,M_opt(TrueStates(t-1),:,k));
            MatchProp = Obs(TrueStates(t))/NumExp;
        end
        MatchPropSum(k) = MatchPropSum(k) + MatchProp;
    end
    NumTransitions(k) = NumTransitions(k) + t_k;
end
% MeanPredMatchEachPatient = MatchPropSum./NumTransitions;
MeanPredMatchOverall = sum(MatchPropSum)/sum(NumTransitions);

end
        
        
    
    

