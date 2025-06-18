function LogLike = log_likelihood_with_BetaMatSelected_FASTEST(SeqCell, X, BetaMatSelected, N,NonNullPositions, MHat)
K = size(X,1); % Number of patients
p = size(X, 2) - 1; % number of covariates (excluding intercept)

if(sum(sum(NonNullPositions)) ~= size(BetaMatSelected,1))
    fprintf('\n');
    fprintf('ERROR: Total number of non-zero locations doesnt match with the the number of rows of BetaMatSelected. \n');
    fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding patientSpecific Trans matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = NonNullPositions; % Indicator matrix                       
J = ones(N+1,N) - I; % complementary indicator matrix

%%%%%%%% Linear projection matrices (patient-specific) %%%%%%%%%%%%%%%%%%%%

G = MHat .* J;    % The empirical estimates                                
GRowSums = sum(G, 2);
HRowBudget = ones(N+1, 1) - GRowSums;

L = zeros(N+1, N, K);
H = zeros(N+1, N, K);
M = zeros(N+1, N, K);

BetaMatSelectedIdx = 0;
for i = 1:(N+1)
    for j = 1:N
        if(I(i,j) == 1)
            BetaMatSelectedIdx = BetaMatSelectedIdx + 1;
            BetaTemp0 = BetaMatSelected(BetaMatSelectedIdx,:)';                          
            BetaTemp = BetaTemp0 / norm(BetaTemp0);
            for k = 1:K
                L(i,j,k) = exp(X(k,:) * BetaTemp);
            end
        end
    end
end

for k = 1:K
    H(:,:,k) = L(:,:,k) .* I;
    TempMat = (H(:,:,k) ./ sum(H(:,:,k), 2)) .* HRowBudget;
    TempMat(isnan(TempMat)) = 0;
    
    M(:,:,k) = G + TempMat;                                                
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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