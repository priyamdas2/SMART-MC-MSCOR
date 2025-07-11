%  This function takes 
%    (i) treatment sequence in (K * n_k) dimensional cell 
%    (ii) patient covariate values,
%    (iii) Beta array (length = [N+1] * N * [p+1]) as inputs
%  It returns
%   (a) (N+1)*N*(p+1) dimensional Patient specific transition matrices
%  =>=>=>  Caution => (1) beta vector is normalized WITHIN the function




function M = Patient_Trans_Mat_robust(SeqCell, X, BetaVec, N)                        % sum(BetaVec ~= 0)/(6)
%%% Extracting no. of treatments, no. of patients, no. of covariates

K = size(X,1); % Number of patients
p = size(X, 2) - 1; % number of covariates (excluding intercept)
CountTol = CountTol_fun_robust(p); % minimum number of count to model it as 
                            % function of covariates

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
%%%  Writing Transition indicator matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 0 => model -> transprob = zero;
%%% 1 => model -> transprob = a function of covariates

InitialIndMat = ones(1, N);
TransIndMat = ones(N, N);

for j = 1:N
    if(IntialCountMat(1,j) <= CountTol)
        InitialIndMat(1,j) = 0;
    end
end

for i = 1:N
    for j = 1:N
        if(TransCountMat(i, j) <= CountTol)
            TransIndMat(i,j) = 0;
        end
    end
end

I = [InitialIndMat; TransIndMat]; % Indicator matrix                       % sum(sum(I))

J = ones(N+1,N) - I; % complementary indicator matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Writing Empirical transition matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MHat = nan(N+1, N); % Empirical transition probabilities

for i = 1:(N+1)
    MHat(i, :) = CHat(i, :) / sum(CHat(i, :));
end

                                                                           % sum(MHat, 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Linear projection matrices (patient-specific) %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = MHat .* J;    % The empirical estimates                                % sum(sum(G == 0))
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
            BetaTemp0 = BetaVec(StartInd:EndInd);                          % = BetaCell{i,j}
            BetaTemp = BetaTemp0 / norm(BetaTemp0);
        end
        for k = 1:K
            if(I(i,j) == 1)
                L(i,j,k) = exp(X(k,:) * BetaTemp);                         % sum(sum((L(:,:,1) ~= 0)))
            end
        end
    end
end

for k = 1:K
    H(:,:,k) = L(:,:,k) .* I;
    TempMat = (H(:,:,k) ./ sum(H(:,:,k), 2)) .* HRowBudget;
    TempMat(isnan(TempMat)) = 0;
    
    M(:,:,k) = G + TempMat;                                                % sum(sum(M(:,:,k), 2))   
end
end



