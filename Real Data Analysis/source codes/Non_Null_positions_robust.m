% It returns the N * N indicator matrix denoting which positions are to be 
% estimated via beta

function [CHat, NonNullPositions] = Non_Null_positions_robust(SeqCell, N, p)
%%% Extracting no. of treatments, no. of patients, no. of covariates
K = max(size(SeqCell)); % Number of patients
CountTol = CountTol_fun_robust(p);

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

NonNullPositions = [InitialIndMat; TransIndMat]; % Indicator matrix
end