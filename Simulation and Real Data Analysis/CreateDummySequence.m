function DummySequences = CreateDummySequence(RandSeed, XDummyPatient, ...
    NumPatientsEachClass, EachSeqLength, SeqCell, X, BetaCell, N)

K = size(X,1); % Number of patients
NumDummyClasses = size(XDummyPatient,1); % Number of Dummy patients
p = size(X, 2) - 1; % number of covariates (excluding intercept)
CountTol = CountTol_fun(p); % minimum number of count to model it as function of covariates


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

L_dummy = zeros(N+1, N, NumDummyClasses);
H_dummy = zeros(N+1, N, NumDummyClasses);
M_dummy = zeros(N+1, N, NumDummyClasses);


BetaVecWhole = Beta_cell2vec(BetaCell);
NonNullPositions = I;
BetaVec = BetaVec_NonNULL_selection(BetaVecWhole, NonNullPositions);

for i = 1:(N+1)
    for j = 1:N
        if(I(i,j) == 1)
            StartInd = (i-1)*N*(p+1)+(j-1)*(p+1) + 1;
            EndInd = (i-1)*N*(p+1)+(j-1)*(p+1) + p + 1;
            BetaTemp0 = BetaVec(StartInd:EndInd);                          
            BetaTemp = BetaTemp0 / norm(BetaTemp0);
        end
        for k = 1:NumDummyClasses
            if(I(i,j) == 1)
                L_dummy(i,j,k) = exp(XDummyPatient(k,:) * BetaTemp);                      
            end
        end
    end
end

for k = 1:NumDummyClasses
    H_dummy(:,:,k) = L_dummy(:,:,k) .* I;
    TempMat = (H_dummy(:,:,k) ./ sum(H_dummy(:,:,k), 2)) .* HRowBudget;
    TempMat(isnan(TempMat)) = 0;
    
    M_dummy(:,:,k) = G + TempMat;                                              
end

rng(RandSeed)

DummySequences = zeros(NumPatientsEachClass, EachSeqLength, NumDummyClasses);

for k = 1:NumDummyClasses
    M_here = M_dummy(:,:,k);
    ISV = M_here(1,:);
    TM = M_here(2:end,:);
    TempSeqMat = zeros(NumPatientsEachClass, EachSeqLength);
    for ii = 1:NumPatientsEachClass
        
        for t = 1:EachSeqLength
            if(t == 1)
                TempSeqMat(ii,t) = find(mnrnd(1, ISV, 1));
            else
                TempSeqMat(ii,t) = find(mnrnd(1, TM(TempSeqMat(ii,t-1),:), 1));
            end
        end
    end
    DummySequences(:,:,k) = TempSeqMat;
end

end


                     
