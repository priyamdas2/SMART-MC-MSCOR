clear all
DoBootstrapNow = 0; % 0 = Bootstrap is performed now; 1 = Bootstrap results are loaded from saved results
%%% MSCOR parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MaxTime = 3600;
MaxRuns = 10;
MaxIter = 5000;
TolFun1 = 10^(-1);
TolFun2 = 10^(-1);
phi = 10^(-20);
lambda = 10^(-6);
RandSeed = 2;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Reading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - Covariates: c(Intercept, "age_at_diag", "disease_duration", "gender","white", "black")
% - Treatments: 
% 1 = 121191 = B-cell depletion (rituximab + ocrelizumab)
% 2 = 214582 = glatiramer acetate
% 3 = NA = Interferon-beta
% 4 = 1373478 = dimethyl fumarate
% 5 = 354770 = natalizumab
% 6 = 1012892 = fingolimod
% 7 = 1310520 = teriflunomide
% 8 = 92299 = cyclophosphamide
% 9 = 7005 = mitoxantrone
% 10 = 117055 = alemtuzumab

%%%  Reading Patient treatment sequence data 
filename = ['Real Data/MS treatment sequences.csv'];
TrtSeqMat = csvread(filename);

%%%  Reading Patient covariate data
filename = ['Real Data/MS covariates.csv'];
XRaw = csvread(filename);

%%% Scaliing X covaraites (except intercept) to (mead,sd) = (0,1)
% While scaling, 'age_at_diag' (X(:,2)) and 'disease_duration' ((X(:,3)) 
% are scaled with (mean = 12, sd = 10) and (mean = 0, sd = 10), Avg. age of
% diagnosis of MS in USA is 34.
% respectively. Other variables are kept unchanged.

% X = [ones(size(XRaw,1), 1), (XRaw(:,2:end) - mean(XRaw(:,2:end)))./ ...
%     std(XRaw(:,2:end))];
age_mu = 12;
age_sigma = 10;
dd_mu = 0;
dd_sigma = 10;

filename = 'Output_age_muSigma_dd_muSigma.csv';
writematrix([age_mu, age_sigma, dd_mu, dd_sigma], filename);

AgeAtDiag_scaled = (XRaw(:,2) - age_mu)./age_sigma;
DiseaseDuration_scaled = (XRaw(:,3) - dd_mu)./dd_sigma;
X = [ones(size(XRaw,1), 1), AgeAtDiag_scaled, DiseaseDuration_scaled,...
    XRaw(:,4:end)];


%%% Extracting no. of treatments, no. of patients, no. of covariates

N = size(unique(TrtSeqMat(:, 2)), 1); % Number of treatments
K = size(unique(TrtSeqMat(:, 1)), 1); % Number of patients
p = size(XRaw, 2) - 1; % number of covariates (excluding intercept)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding Bootstrap Standard errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumBootRep = 1000; % Dont change; saved result is available only for value 1000

if(DoBootstrapNow == 1)
    BetaALLBootStdErr = BetaALL_BootStdErr(NumBootRep, X, TrtSeqMat, N);
    filename = ['BetaALL_REALDATA_BootStdErr_rep_',num2str(NumBootRep),'.mat'];
    save(filename, 'BetaALLBootStdErr');
else
    filename = ['BetaALL_REALDATA_BootStdErr_rep_',num2str(NumBootRep),'.mat'];
    TempMat = load(filename);
    BetaALLBootStdErr = TempMat.BetaALLBootStdErr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Trying out different Important functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  Converting trt-seq data from matrix to cell format
SeqCell = Trt_Seq_Mat2Cell(TrtSeqMat, N);

%%% Identifying positions to be estimated by beta
[CHat, NonNullPositions] = Non_Null_positions(SeqCell, N, p);
MHat = CHat2MHat(CHat);

filename = 'Output_CHat.csv';
writematrix(CHat, filename);

filename = 'Output_NonNullPositions.csv';
writematrix(NonNullPositions, filename);

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
    1, 2, 1, 0);


%%% Printing Initial and Final 
Initial_LogLikelihood = -objFun(x0MatSelected);
Final_LogLikelihood = -objFun(x_opt_par);

fprintf('Log-likelihood at starting point: %.2f \n',...
    Initial_LogLikelihood); 
fprintf('Log-likelihood at solution point: %.2f \n',...
    Final_LogLikelihood); 

BetaShortVecUnscaled_opt = BetaMatSelected2ShortVec(x_opt_par); 

%%%% MSCORparallel solution is normalized already, normalizing once more 
% to make norm exactly 0 with high precision %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BetaShortVecNormalized_opt = BetaShortVec_normalizer(BetaShortVecUnscaled_opt, p);

BetaShortVec_opt = BetaShortVecNormalized_opt;
BetaVec_opt = BetaShortVec2Vec(BetaShortVec_opt,NonNullPositions);

filename = 'Output_BetaVec_opt.csv';
writematrix(BetaVec_opt, filename);

BetaCell_opt = Beta_vec2cell(BetaVec_opt, N);
M_opt = Patient_Trans_Mat(SeqCell, X, BetaVec_opt, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Comparison with NAIVE MC model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Should NOT expect much improvement

MHat = MHatEstimate(SeqCell, X, N);
M_naive = zeros(N+1,N,K);
for k = 1:K
    M_naive(:,:,k) = MHat;
end


Naive_MeanPredMatchOverall = MeanPredMatch(M_naive, SeqCell, 56, 1000);
Proposed_MeanPredMatchOverall = MeanPredMatch(M_opt, SeqCell, 56, 1000);
[Naive_MeanPredMatchOverall, Proposed_MeanPredMatchOverall]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Patient specific transition matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Covariates HERE: c(Intercept, "age_at_diag", "disease_duration", 
%    "gender","white", "black")
% find((XRaw(:, 2) > 20) .* (XRaw(:, 2) < 30) .* (XRaw(:, 5) == 1))' # id 17
% find((XRaw(:, 2) > 20) .* (XRaw(:, 2) < 25) .* (XRaw(:, 6) == 1))' # id 412
% find((XRaw(:, 2) > 40) .* (XRaw(:, 2) < 50) .* (max(XRaw(:, 5), XRaw(:,6)) == 0))' # id 674
% find((XRaw(:, 2) > 60) .* (XRaw(:, 2) < 90))'  # id 796

MFinalCell = M_opt;
PatientNumbers = [17, 412, 674, 796];
Plot_PatientTransitionProb(MFinalCell, XRaw, PatientNumbers)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta Coefficients  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BetaIndex : 2 = age_at_diag, 3 = disease_duration, 4 = gender",
%             5 = white, 6 = black.

BetaCellFinal = BetaCell_opt;
for BetaIndex = 2:6
    Plot_Beta(BetaIndex,BetaCellFinal,NonNullPositions)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Most frequent treatment transition analyses  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TopHowMany = 15;
[CoeffSummary, SeSummary] = TopHowManyTransFormula(TopHowMany, SeqCell, BetaCellFinal,...
    BetaALLBootStdErr, N);

filename = 'Output_CoeffSummary.csv';
writematrix(CoeffSummary, filename);

filename = 'Output_SeSummary.csv';
writematrix(SeSummary, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Patient-specific example treatment sequence %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumDummyPatients = 12;
XDummyPatient_Raw = zeros(NumDummyPatients, (p+1));
% mean(XRaw(:,3)) = 15.4063 % mean disease duration
XDummyPatient_Raw(1,:) = [1 30 15.4063 1 1 0]; % (Age = 30, F, White)
XDummyPatient_Raw(2,:) = [1 30 15.4063 0 1 0]; % (Age = 30, M, White)
XDummyPatient_Raw(3,:) = [1 30 15.4063 1 0 1]; % (Age = 30, F, Black)
XDummyPatient_Raw(4,:) = [1 30 15.4063 0 0 1]; % (Age = 30, M, Black)
XDummyPatient_Raw(5,:) = [1 30 15.4063 1 0 0]; % (Age = 30, F, Others)
XDummyPatient_Raw(6,:) = [1 30 15.4063 0 0 0]; % (Age = 30, M, Others)
XDummyPatient_Raw(7,:) = [1 60 15.4063 1 1 0]; % (Age = 60, F, White)
XDummyPatient_Raw(8,:) = [1 60 15.4063 0 1 0]; % (Age = 60, M, White)
XDummyPatient_Raw(9,:) = [1 60 15.4063 1 0 1]; % (Age = 60, F, Black)
XDummyPatient_Raw(10,:) = [1 60 15.4063 0 0 1]; % (Age = 60, M, Black)
XDummyPatient_Raw(11,:) = [1 60 15.4063 1 0 0]; % (Age = 60, F, Others)
XDummyPatient_Raw(12,:) = [1 60 15.4063 0 0 0]; % (Age = 60, M, Others)

%%% Scaling 
AgeAtDiag_scaled = (XDummyPatient_Raw(:,2) - age_mu)./age_sigma;
DiseaseDuration_scaled = (XDummyPatient_Raw(:,3) - dd_mu)./dd_sigma;
XDummyPatient = [ones(size(XDummyPatient_Raw,1), 1), AgeAtDiag_scaled, DiseaseDuration_scaled,...
    XDummyPatient_Raw(:,4:end)];

%%% Generating treatment sequences

RandSeed = 1;
EachSeqLength = 20;
NumPatientsEachClass = 50;
DummySequences = CreateDummySequence(RandSeed, XDummyPatient, ...
    NumPatientsEachClass, EachSeqLength, SeqCell, X, BetaCell_opt, N);

%%% 
Plot_DummySequences_20Treatments(DummySequences,XDummyPatient_Raw)