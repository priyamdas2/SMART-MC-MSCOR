function LogLike = log_likelihood_NaiveMethod(SeqCell, X, BetaVec, N)
K = size(X,1); % Number of patients
p = size(X, 2) - 1; % number of covariates (excluding intercept)
M = Patient_Trans_Mat_NaiveMethod(X, BetaVec, N);


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