function CHat = TransitionCountMat(SeqCell, N)                        % sum(BetaVec ~= 0)/(6)
%%% Extracting no. of treatments, no. of patients, no. of covariates

K = size(SeqCell, 2); % Number of patients

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
end