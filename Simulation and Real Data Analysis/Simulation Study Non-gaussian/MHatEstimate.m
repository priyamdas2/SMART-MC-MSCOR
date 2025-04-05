function MHat = MHatEstimate(SeqCell, X, N)                        % sum(BetaVec ~= 0)/(6)
%%% Extracting no. of patients, no. of covariates
K = size(X,1); % Number of patients
p = size(X, 2) - 1; % number of covariates (excluding intercept)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Finding transition frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IntialCountMat = zeros(1, N);
TransCountMat = zeros(N, N);
for i = 1:K
    
    for j = 1:(length(SeqCell{i})-1)
        aa = SeqCell{i}(j);
        bb = SeqCell{i}(j+1);
        if(j == 1)
            IntialCountMat(1,aa) = IntialCountMat(1,aa)+1;
        end
        TransCountMat(aa,bb) = TransCountMat(aa,bb) + 1;
    end
end
CHat = [IntialCountMat; TransCountMat];  % Empirical count matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Writing Empirical transition matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MHat = nan(N+1, N); % Empirical transition probabilities

for i = 1:(N+1)
    MHat(i, :) = CHat(i, :) / sum(CHat(i, :));
end