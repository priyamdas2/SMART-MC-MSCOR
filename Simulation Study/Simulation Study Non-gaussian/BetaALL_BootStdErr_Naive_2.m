function BetaALLBootStdErr = BetaALL_BootStdErr_Naive (NumBootRep, X, TrtSeqMat, N)

K = size(X ,1); % Number of patients
p = size(X, 2) - 1; % number of covariates (excluding intercept)

BetaALLBoot = cell(N+1,N,p+1);
SeqCell = Trt_Seq_Mat2Cell(TrtSeqMat, N);

time_spents = nan(NumBootRep,1);
for iter = 1:NumBootRep
    rng(iter)
    tic;
    fprintf('Performing Bootstrap Iteration: %d  \n',iter);
    RandSeed = iter;
    population = 1:K;  % Example population
    num_samples = K;
    sample_indices = randi(length(population), 1, num_samples);  % Random indices with replacement
    samples = population(sample_indices);  % Extract samples using the indices
    X_Boot = X(samples,:);
    
    SeqCell_Boot = cell(1,K);
    for k = 1:K
        TrtSeqLength = length(SeqCell{samples(k)});
        SeqCell_Boot{k}(1:TrtSeqLength) = SeqCell{samples(k)}(1:TrtSeqLength);
    end
        
    objFunc_Boot = @(z) -log_likelihood_NaiveMethod(SeqCell, X_Boot, z, N);
    
    x0Mat_Boot = Generate_Rand_Initial_BetaMat_2(RandSeed, N, p);    
    x0_Boot = BetaMat2Vec(x0Mat_Boot);
    d = length(x0_Boot);  
    lb = -1*ones(d,1);  % Lower bounds
    ub = 1*ones(d,1);  % Upper bounds
    
    options = optimoptions('fmincon','MaxFunctionEvaluations',10000,...
    'Display', 'off');
    [x_opt_Boot, fval] = fmincon(objFunc_Boot, x0_Boot, [], [], [], [], lb, ub,[], options);
    BetaVec_opt_Boot = x_opt_Boot;
    BetaCell_opt_Boot = Beta_vec2cell(BetaVec_opt_Boot, N);
    
    %%% Normalizing each block to have L2 norm 1
    BetaCell_opt_Boot_scaled = BetaCell_opt_Boot;  % Initialize scaled version
    
    [nRows, nCols] = size(BetaCell_opt_Boot);  % Get cell array size
    
    for i = 1:nRows
        for j = 1:nCols
            vec = BetaCell_opt_Boot{i, j};  % Extract vector
            norm_val = norm(vec, 2);  % Compute L2 norm
            
            if norm_val > 0
                BetaCell_opt_Boot_scaled{i, j} = vec / norm_val;  % Normalize
            else
                BetaCell_opt_Boot_scaled{i, j} = vec;  % Keep zero vector unchanged
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:(N+1)
        for j = 1:N
            for d = 1:(p+1)
                if(iter == 1)
                    BetaALLBoot{i,j,d} = BetaCell_opt_Boot_scaled{i,j}(d);
                else
                    BetaALLBoot{i,j,d} = [BetaALLBoot{i,j,d}; BetaCell_opt_Boot_scaled{i,j}(d)];
                end
            end
        end
    end
    time_spents(iter) = toc;
    fprintf('\n');
    fprintf('=> Time spent on last Bootstrap iteration: %.2f seconds.  \n',time_spents(iter));
    fprintf('\n');
    EstTimeNeeded = (NumBootRep - iter)*mean(time_spents(1:iter))/60;
    fprintf('=> Estimated additional time to complete: %.2f minutes.  \n',EstTimeNeeded);
    fprintf('\n');
end


BetaALLBootStdErr = nan(N+1,N,p+1);
for i = 1:(N+1)
    for j = 1:N
        for d = 1:(p+1)
            Vec = BetaALLBoot{i,j,d};
            nonNaN_idx = ~isnan(Vec);
            Vec_nonNaN = Vec(nonNaN_idx);
            Vec_NZ = Vec_nonNaN(Vec_nonNaN ~=0);
            BetaALLBootStdErr(i,j,d) = std(Vec_NZ)/sqrt(length(Vec_NZ));
        end
    end
end

end
