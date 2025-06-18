clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Data generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% MSCOR parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MaxTime = 3600;
MaxRuns = 10;
MaxIter = 5000;
TolFun1 = 10^(-1);
TolFun2 = 10^(-1);
phi = 10^(-20);
lambda = 10^(-6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumExperiment = 10;
ALL_TrtSeqLengthEach = [20,40,60];
ALL_NumSubjects = [1000,2000,3000];
percentage = 10;
n1 = length(ALL_TrtSeqLengthEach);
n2 = length(ALL_NumSubjects);

for repIdx = 1:NumExperiment
    RandSeed = repIdx;
    ALL_MADs = nan(n1,n2);
    for rowIdx = 1:n1
        for colIdx = 1:n2
            
            TrtSeqLengthEach = ALL_TrtSeqLengthEach(rowIdx);
            NumSubjects = ALL_NumSubjects(colIdx);
            
            fprintf('\n');
            fprintf('=> Now performing exp. no.: %d, scenario (%d,%d): Seq length %d, number of subjects %d. \n',repIdx, rowIdx, colIdx, TrtSeqLengthEach,NumSubjects);
            fprintf('\n');
            
            NumContVars = 5;
            NumIndicVars = 0;
            InterceptIndic = 1;
            mu = 0;            % Used in Beta generation
            sigma = 10;        % Used in Beta generation
            BernProb = 0.5;    % Used if dist = 1, and, NumIndicVars > 0
            dist = 2;          % X generated from: 1 = uniform(0,1), 2 = Normal(0,1)
            
            p = NumContVars+NumIndicVars;
            N = 10;             % use >= 5
            TrueNullPercentage = 67;  % This many non-diagonal transition probabilities are made zero
            
            NonNullPositionsTrue = Simulation_NonNullPositionMatrixGeneratorCORRECTED(RandSeed,TrueNullPercentage,N);
            
            X = Simulation_generate_X(RandSeed,NumSubjects,NumContVars,NumIndicVars,...
                BernProb,InterceptIndic,dist);
            
            BetaCellTrue = BetaCellRandGenerator(RandSeed, N, p, mu, sigma);
            SubjectTransMats = SubjectSpecificTransMat(X,BetaCellTrue,NonNullPositionsTrue);
            TrtSeqMat = Simulation_generate_TrtSeqMat(RandSeed,SubjectTransMats,...
                TrtSeqLengthEach);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%% SMART-MC estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            SeqCell = Trt_Seq_Mat2Cell(TrtSeqMat, N);
            [CHat, NonNullPositions] = Non_Null_positions(SeqCell, N, p);
            MHat = CHat2MHat(CHat);
            
            %%% Random initial starting point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x0Mat = Generate_Rand_Initial_BetaMat(RandSeed, N, p);
            x0MatSelected = BetaMat2Select(x0Mat,N, NonNullPositions);
            
            %%% Optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % log_likelihood_with_BetaMatSelected_FASTEST(SeqCell, X, x0MatSelected, N,NonNullPositions, MHat)
            
            objFun = @(zMat) -log_likelihood_with_BetaMatSelected_FASTEST(SeqCell,...
                X, zMat, N,NonNullPositions, MHat);
            
            %[x_opt, fval, comp_time] = MSCOR(fun, x0MatSelected, MaxTime, MaxRuns, MaxIter, sInitial, rho, TolFun1, TolFun2, phi,lambda, DisplayUpdate, DisplayEvery,PrintStepSize, PrintSolution)
            [x_opt_par, fval_par, comp_time_par] = MSCORparallel(objFun, x0MatSelected, MaxTime, MaxRuns, MaxIter, 1, 2, TolFun1, TolFun2, phi, lambda, 1, 2, 1, 0);
            BetaShortVecUnscaled_opt = BetaMatSelected2ShortVec(x_opt_par);
            
            %%%% MSCORparallel solution is normalized already, normalizing once more
            % for perfection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BetaShortVecNormalized_opt = BetaShortVec_normalizer(BetaShortVecUnscaled_opt, p);
            BetaShortVec_opt = BetaShortVecNormalized_opt;
            BetaVec_opt = BetaShortVec2Vec(BetaShortVec_opt,NonNullPositions);
            BetaCell_opt = Beta_vec2cell(BetaVec_opt, N);
            BetaCellEst = BetaCell_opt;
            
            TopHowManyInitials = round((percentage/100)*N);
            TopHowManyTrans = round((percentage/100)*(N^2));
            % TopHowManyInitials = 2;
            % TopHowManyTrans = 5;
            [CoeffsEst, CoeffsTrue, DiffMatTrans, RMSE, MAD]=Simulation_TopHowManyBetaComparison(TopHowManyInitials,...
                TopHowManyTrans, SeqCell,BetaCellEst,BetaCellTrue,InterceptIndic,N);
            ALL_MADs(rowIdx,colIdx)= MAD;
        end
    end
    ALL_MADs_2 = [ALL_TrtSeqLengthEach',ALL_MADs];
    ALL_MADs_3 = [[nan, ALL_NumSubjects];ALL_MADs_2];
    
    
    filename = ['MAD outputs/Output_Simulation_MAD_sequence_RandSeed_',num2str(repIdx),'.csv'];
    writematrix(round(ALL_MADs_3,3), filename);
end



Temp = cell(NumExperiment,1);
for repIdx = 1:NumExperiment
    filename = ['MAD outputs/Output_Simulation_MAD_sequence_RandSeed_',num2str(repIdx),'.csv'];
    Temp{repIdx} = readmatrix(filename);
end

SummaryMean = nan(n1,n2);
SummarySe = nan(n1,n2);
for ii = 1:n1
    for jj = 1:n2
        for repIdx = 1:NumExperiment
            if(repIdx == 1)
                valueArray = Temp{repIdx}(ii+1,jj+1);
            else
                valueArray = [valueArray;Temp{repIdx}(ii+1,jj+1)];
            end
        end
        SummaryMean(ii,jj) = mean(valueArray);
        SummarySe(ii,jj) = std(valueArray)/sqrt(length(valueArray));
    end
end

filename = ['Output_Simulation_MAD_Mean.csv'];
writematrix(round(SummaryMean,4), filename);

filename = ['Output_Simulation_MAD_Se.csv'];
writematrix(round(SummarySe,4), filename);

