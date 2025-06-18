function LogLike = log_likelihood_with_BetaMat_FASTER(SeqCell, X, BetaMat, N,NonNullPositions, MHat)
K = size(X,1); % Number of patients
p = size(X, 2) - 1; % number of covariates (excluding intercept)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Converting BetaMat to BetaVec %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumRows = size(BetaMat,1);
for i = 1:NumRows
    if(i == 1)
        BetaVec = BetaMat(i,:)';
    else
        BetaVec = [BetaVec; BetaMat(i,:)'];
    end
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


for i = 1:(N+1)
    for j = 1:N
        if(I(i,j) == 1)
            StartInd = (i-1)*N*(p+1)+(j-1)*(p+1) + 1;
            EndInd = (i-1)*N*(p+1)+(j-1)*(p+1) + p + 1;
            BetaTemp0 = BetaVec(StartInd:EndInd);                          
            BetaTemp = BetaTemp0 / norm(BetaTemp0);
        end
        for k = 1:K
            if(I(i,j) == 1)
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