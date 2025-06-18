clear all

%%% MSCOR parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MaxTime = 3600;
MaxRuns = 10;
MaxIter = 5000;
TolFun1 = 10^(-1);
TolFun2 = 10^(-1);
phi = 10^(-20);
lambda = 10^(-6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Data generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DoBootstrapNow = 1; % 1 = Bootstrap is performed now; 0 = Bootstrap s.e. are loaded from saved results
RandSeed = 2; % Don't change, 'NonNullPositions' depends on it. And hence calculated bootstrap s.e.
NumSubjects = 1000;
NumContVars = 5;
NumIndicVars = 0;
InterceptIndic = 1;
mu = 0;            % Used in Beta generation
sigma = 10;        % Used in Beta generation
BernProb = 0.5;    % Used if dist = 1, and, NumIndicVars > 0
dist = 2;          % X generated from: 1 = uniform(-10,10), 2 = Normal(0,1)
TrtSeqLengthEach = 20;
p = NumContVars+NumIndicVars;
N = 10;             % use >= 5 
TrueNullPercentage = 67;  % This many non-diagonal transition probabilities are made zero

NonNullPositionsTrue = Simulation_NonNullPositionMatrixGeneratorCORRECTED(RandSeed,TrueNullPercentage,N);

X = Simulation_generate_X(RandSeed,NumSubjects,NumContVars,NumIndicVars,...
    BernProb,InterceptIndic,dist);

BetaCellTrue = BetaCellRandGenerator(RandSeed, N, p, mu, sigma);
save('BetaCellTrue.mat', 'BetaCellTrue')

SubjectTransMats = SubjectSpecificTransMat(X,BetaCellTrue,NonNullPositionsTrue);
TrtSeqMat = Simulation_generate_TrtSeqMat(RandSeed,SubjectTransMats,...
    TrtSeqLengthEach);

writematrix(X, 'Dummy_X.csv');
writematrix(TrtSeqMat, 'Dummy_sequences.csv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding Bootstrap Standard errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumBootRep = 100; % Dont change; saved result is available only for value 100

if(DoBootstrapNow == 1)
    BetaALLBootStdErr = BetaALL_BootStdErr(NumBootRep, X, TrtSeqMat, N);
    filename = ['BetaALL_SIMULATION_BootStdErr_rep_',num2str(NumBootRep),'.mat'];
    save(filename, 'BetaALLBootStdErr');
else
    filename = ['BetaALL_SIMULATION_BootStdErr_rep_',num2str(NumBootRep),'.mat'];
    TempMat = load(filename);
    BetaALLBootStdErr = TempMat.BetaALLBootStdErr;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SMART-MC estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SeqCell = Trt_Seq_Mat2Cell(TrtSeqMat, N);
[CHat, NonNullPositions] = Non_Null_positions(SeqCell, N, p);
MHat = CHat2MHat(CHat);

%%% Random initial starting point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0Mat = Generate_Rand_Initial_BetaMat(RandSeed, N, p);
x0MatSelected = BetaMat2Select(x0Mat,N, NonNullPositions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Random initial starting point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0Mat = Generate_Rand_Initial_BetaMat(RandSeed, N, p);
x0MatSelected = BetaMat2Select(x0Mat,N, NonNullPositions);

log_likelihood_with_BetaMatSelected_FASTEST(SeqCell, X, x0MatSelected,N,...
    NonNullPositions, MHat)

objFun = @(zMat) -log_likelihood_with_BetaMatSelected_FASTEST(SeqCell,...
    X, zMat, N,NonNullPositions, MHat);

%[x_opt, fval, comp_time] = MSCOR(fun, x0MatSelected, MaxTime, MaxRuns, MaxIter, sInitial, rho, TolFun1, TolFun2, phi,lambda, DisplayUpdate, DisplayEvery,PrintStepSize, PrintSolution)    
[x_opt_par, fval_par, comp_time_par] = MSCORparallel(objFun, x0MatSelected,...
    MaxTime, MaxRuns, MaxIter, 1, 2, TolFun1, TolFun2, phi, lambda,...
    0, 2, 1, 0);
BetaShortVecUnscaled_opt = BetaMatSelected2ShortVec(x_opt_par); 

%%%% MSCORparallel solution is normalized already, normalizing once more 
% to make norm exactly 0 with high precision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BetaShortVecNormalized_opt = BetaShortVec_normalizer(BetaShortVecUnscaled_opt, p);

BetaShortVec_opt = BetaShortVecNormalized_opt;
BetaVec_opt = BetaShortVec2Vec(BetaShortVec_opt,NonNullPositions);
BetaCell_opt = Beta_vec2cell(BetaVec_opt, N);
BetaCellEst = BetaCell_opt;


NonNullPositionsTrue - NonNullPositions;
percentage = 10;
TopHowManyInitials = round((percentage/100)*N);
TopHowManyTrans = round((percentage/100)*(N^2));
[CoeffsEst, CoeffsTrue, DiffMatTrans, RMSE, MAD]=Simulation_TopHowManyBetaComparison(TopHowManyInitials,...
    TopHowManyTrans, SeqCell,BetaCellEst,BetaCellTrue,InterceptIndic,N);
    
CoeffsEst(:,5:end)
CoeffsTrue(:,5:end)
DiffMatTrans
MAD

filename = 'Output_Simulation_CoeffsEst.csv';
writematrix(round(CoeffsEst,2), filename);

filename = 'Output_Simulation_CoeffsTrue.csv';
writematrix(round(CoeffsTrue,2), filename);

% CoeffsEst = readmatrix('Output_Simulation_CoeffsEst.csv')
% CoeffsTrue = readmatrix('Output_Simulation_CoeffsTrue.csv')
% CoeffsDiffMAD = mean(abs(CoeffsEst(:,5:10) - CoeffsTrue(:,5:10)), 2);
% filename = 'Output_Simulation_CoeffsTrueEstMAD.csv';
% writematrix(round(CoeffsDiffMAD,2), filename);


TopTransitions = CoeffsEst(:,1:2);
Se = zeros(TopHowManyInitials+TopHowManyTrans, NumContVars + 1);

for ii = 1:(TopHowManyInitials+TopHowManyTrans)
    if(ii <= TopHowManyInitials)
        coeffs_here = BetaALLBootStdErr(1,TopTransitions(ii,1), 1);
        for coeffIdx = 2:(NumContVars + 1)
            coeffs_here = [coeffs_here, BetaALLBootStdErr(1,TopTransitions(1,1), coeffIdx)];
        end
        Se(ii,:) = coeffs_here;
    else
        coeffs_here = BetaALLBootStdErr(TopTransitions(ii,1) + 1,TopTransitions(ii,2), coeffIdx);
        for coeffIdx = 2:(NumContVars + 1)
            coeffs_here = [coeffs_here, BetaALLBootStdErr(TopTransitions(ii,1) + 1,TopTransitions(ii,2), coeffIdx)];
        end
        Se(ii,:) = coeffs_here;
    end
end

filename = 'Output_Simulation_StdErr.csv';
writematrix(round(Se,3), filename);


