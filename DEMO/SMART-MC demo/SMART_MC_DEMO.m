clear all
% Takes ~ 218 seconds on Windows 10 Enterprise operating system with
% 32 GB RAM and the following processor specifications: 12th Gen Intel(R)
% Core(TM) i7-12700, 2100 MHz, 12 Cores, and 20 Logical Processors.
RandSeed = 1;
%% MSCOR parameters

ActivateMSCORParallel = 1; % 1 = parallel, 0 = NOT parallel
MaxTime = 3600;
MaxRuns = 10;
MaxIter = 5000;
sInitial = 1;
rho = 2;
TolFun1 = 10^(-0);
TolFun2 = 10^(-0);
phi = 10^(-4);
lambda = 10^(-20);
DisplayUpdate = 1;
DisplayEvery = 2;  % displays update every this many seconds
PrintStepSize = 1; % displays current step-size (0/1)
PrintSolution = 0; % prints final solution (0/1)


%% Reading data

filename = ['Synthetic Data/Dummy_sequences.csv'];
TrtSeqMat = csvread(filename);

filename = ['Synthetic Data/Dummy_X.csv'];      % This X is scaled.
X = csvread(filename);

%% Extracting no. of treatments, no. of patients, no. of covariates

N = size(unique(TrtSeqMat(:, 2)), 1); % Number of treatments
K = size(unique(TrtSeqMat(:, 1)), 1); % Number of patients
p = size(X, 2) - 1; % number of covariates (excluding intercept)

%% Pre-analysis processing

%%%  Converting trt-seq data from matrix to cell format
SeqCell = Trt_Seq_Mat2Cell(TrtSeqMat, N);

%%% Identifying positions to be estimated by beta
[CHat, NonNullPositions] = Non_Null_positions(SeqCell, N, p);
MHat = CHat2MHat(CHat);

%% MSCOR

%%% Random initial starting point 

x0Mat = Generate_Rand_Initial_BetaMat(RandSeed, N, p);
x0MatSelected = BetaMat2Select(x0Mat,N, NonNullPositions);

objFun = @(zMat) -log_likelihood_with_BetaMatSelected_FASTEST(SeqCell,...
    X, zMat, N,NonNullPositions, MHat);

if(ActivateMSCORParallel == 1)
    [x_opt_par, fval_par, comp_time_par] = MSCORparallel(objFun, x0MatSelected,...
        MaxTime, MaxRuns, MaxIter, sInitial, rho, TolFun1, TolFun2, phi, lambda,...
        DisplayUpdate, DisplayEvery, PrintStepSize, PrintSolution);
else
    [x_opt_par, fval_par, comp_time_par] = MSCOR(objFun, x0MatSelected,...
        MaxTime, MaxRuns, MaxIter, sInitial, rho, TolFun1, TolFun2, phi, lambda,...
        DisplayUpdate, DisplayEvery, PrintStepSize, PrintSolution);
end


%%% Printing Initial and Final log-likelihoods
Initial_LogLikelihood = -objFun(x0MatSelected);
Final_LogLikelihood = -objFun(x_opt_par);

fprintf('Log-likelihood at starting point: %.2f \n',...
    Initial_LogLikelihood);
fprintf('Log-likelihood at solution point: %.2f \n',...
    Final_LogLikelihood);

BetaShortVecUnscaled_opt = BetaMatSelected2ShortVec(x_opt_par);

%%%% MSCORparallel solution is normalized already, normalizing once more
% to make norm exactly 1 with high precision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BetaShortVecNormalized_opt = BetaShortVec_normalizer(BetaShortVecUnscaled_opt, p);

BetaShortVec_opt = BetaShortVecNormalized_opt;
BetaVec_opt = BetaShortVec2Vec(BetaShortVec_opt,NonNullPositions);

BetaCell_opt = Beta_vec2cell(BetaVec_opt, N);          % Estimated Coefficients
M_opt = Patient_Trans_Mat(SeqCell, X, BetaVec_opt, N); % patient-specific transition matrices
