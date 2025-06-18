%  This function takes 
%    (i) treatment sequence in (K * n_k) dimensional cell 
%    (ii) patient covariate values,
%    (iii) Beta array (length = [N+1] * N * [p+1]) as inputs
%  It returns
%   (a) (N+1)*N*(p+1) dimensional Patient specific transition matrices
%  =>=>=>  Caution => (1) beta vector is normalized WITHIN the function




function M = Patient_Trans_Mat_NaiveMethod(X, BetaVec, N)                        % sum(BetaVec ~= 0)/(6)
%%% Extracting no. of treatments, no. of patients, no. of covariates

K = size(X,1); % Number of patients
p = size(X, 2) - 1; % number of covariates (excluding intercept)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Linear projection matrices (patient-specific) %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = zeros(N+1, N, K);
H = zeros(N+1, N, K);
M = zeros(N+1, N, K);


for i = 1:(N+1)
    for j = 1:N
            StartInd = (i-1)*N*(p+1)+(j-1)*(p+1) + 1;
            EndInd = (i-1)*N*(p+1)+(j-1)*(p+1) + p + 1;
            BetaTemp0 = BetaVec(StartInd:EndInd);                          % = BetaCell{i,j}
            BetaTemp = BetaTemp0 / norm(BetaTemp0);
        for k = 1:K
                L(i,j,k) = exp(X(k,:) * BetaTemp);                         % sum(sum((L(:,:,1) ~= 0)))
        end
    end
end

for k = 1:K
    H(:,:,k) = L(:,:,k);
    TempMat = (H(:,:,k) ./ sum(H(:,:,k), 2));
    %TempMat(isnan(TempMat)) = 0;
    M(:,:,k) = TempMat;                                                % sum(sum(M(:,:,k), 2))   
end
end



