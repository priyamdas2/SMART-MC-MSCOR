clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Data generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RandSeed = 1;
RandSeed2 = 2;
MaxTime = 3600*4;
MaxRuns = 10;
MaxIter = 5000;
study = 3;  % Change 1/2/3
if(study == 1)
    TrtSeqLengthEach = 10; NumSubjects = 1000;
else
    if(study == 2)
        TrtSeqLengthEach = 20; NumSubjects = 1000;
    else
        if(study == 3)
            TrtSeqLengthEach = 20; NumSubjects = 2000;
        end
    end
end

ALL_NumContVars = [3,5,8];
ALL_N = [6,9,12]; % >=6 it TrueNullPercentage = 67
n1 = length(ALL_NumContVars);
n2 = length(ALL_N);
Time_summary = nan(n1*n2,7);

for index = 1:(n1*n2)%[1,2,3,4,5,6,7,8,9]
    rowIdx = ceil(index/n2);
    
    if(mod(index,n2) == 0)
        colIdx = n2;
    else
        colIdx = mod(index,n2);
    end
    NumContVars = ALL_NumContVars(rowIdx);
    N = ALL_N(colIdx);
    
    fprintf('\n');
    fprintf('=> Now performing scenario (%d,%d): Seq length %d, number of subjects %d. \n',rowIdx, colIdx, TrtSeqLengthEach,NumSubjects);
    fprintf('\n');
    
    NumIndicVars = 0;
    InterceptIndic = 1;
    mu = 0;            % Used in Beta generation ~ N(mu, sigma)
    sigma = 10;        % Used in Beta generation
    BernProb = 0.5;    % Used if dist = 1, and, NumIndicVars > 0
    dist = 2;          % X generated from: 1 = uniform(0,1), 2 = Normal(0,1)
    
    p = NumContVars+NumIndicVars; 
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
    
    log_likelihood_with_BetaMatSelected_FASTEST(SeqCell, X, x0MatSelected, N,NonNullPositions, MHat)

    objFun = @(zMat) -log_likelihood_with_BetaMatSelected_FASTEST(SeqCell,...
        X, zMat, N,NonNullPositions, MHat);
    
    % [x_opt, fval, comp_time] = MSCOR(fun, x0MatSelected, MaxTime, MaxRuns, MaxIter, sInitial, rho, TolFun1, TolFun2, phi,lambda, DisplayUpdate, DisplayEvery,PrintStepSize, PrintSolution)
    [x_opt, fval, comp_time] = MSCOR(objFun, x0MatSelected, MaxTime, MaxRuns, MaxIter, 1, 2, 10^(-0), 10^(-0), 10^(-4), 10^(-20), 1, 2, 1, 0);
    [x_opt_par, fval_par, comp_time_par] = MSCORparallel(objFun, x0MatSelected, MaxTime, MaxRuns, MaxIter, 1, 2, 10^(-0), 10^(-0), 10^(-4), 10^(-20), 1, 2, 1, 0);
    
    filename = 'Output_Simulation_MSCOR_time_comparison.csv';
    Time_summary_output2 = readmatrix(filename);
    Time_summary_output2(index,:) = [NumContVars, N,MaxRuns,MaxIter, comp_time,comp_time_par];
    writematrix(round(Time_summary_output2,3), filename);
    
    Time_summary(index,:) = [NumContVars, N,MaxRuns,MaxIter, comp_time,...
        comp_time_par, round(comp_time./comp_time_par,1)];
end


filename = ['Output_Simulation_MSCOR_time_comparison_CaseNo_',num2str(study),'.csv'];
writematrix(round(Time_summary), filename);


