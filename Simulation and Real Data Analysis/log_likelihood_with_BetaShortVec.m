% 'BetaShortVec' is normalized WITHIN this function while finding the
% Likelihood value
function LogLike = log_likelihood_with_BetaShortVec(SeqCell, X, BetaShortVec, NonNullPositions, N)
BetaVec = BetaShortVec2Vec(BetaShortVec,NonNullPositions);
K = size(X,1); % Number of patients
M = Patient_Trans_Mat(SeqCell, X, BetaVec, N);


LogLike = 0;

for k = 1:K
    InitialProb = M(1,SeqCell{k}(1),k);
    TransMat = M(2:end,:,k);
    TransProbIndices = N * (SeqCell{k}(2:end) - 1) + SeqCell{k}(1:(end - 1));
    Like_k = InitialProb*prod(TransMat(TransProbIndices));
    if(Like_k == 0)
        Like_k =  10^(-323);
    end
    LogLike = LogLike + log(Like_k);
end
end