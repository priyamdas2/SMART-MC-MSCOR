clearvars; clc; close all;
addpath('./Real Data/');
addpath('./source codes/');
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
% 1 = B-cell depletion (rituximab + ocrelizumab)  (original code = 1)
% 2 = glatiramer acetate                          (original code = 2)
% 3 = Interferon-beta                             (original code = 3)
% 4 = dimethyl fumarate                           (original code = 4)
% 5 = natalizumab                                 (original code = 5)
% 6 = S1P modulators (fingolimod + teriflunomide) (original code = 6,7)
% 7 = Aggressive / Legacy therapies (cyclophosphamide, mitoxantrone, alemtuzumab) (original code = 8, 9, 10)


%%%  Reading Patient treatment sequence data 
filename = ['MS_treatment_sequences_collapsed.csv'];
TrtSeqMat = csvread(filename);

%%%  Reading Patient covariate data
filename = ['MS covariates.csv'];
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
    BetaALLBootStdErr = BetaALL_BootStdErr_robust(NumBootRep, X, TrtSeqMat, N);
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
[CHat, NonNullPositions] = Non_Null_positions_robust(SeqCell, N, p);
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

tic;
%[x_opt, fval, comp_time] = MSCOR(fun, x0MatSelected, MaxTime, MaxRuns, MaxIter, sInitial, rho, TolFun1, TolFun2, phi,lambda, DisplayUpdate, DisplayEvery,PrintStepSize, PrintSolution)    
[x_opt_par, fval_par, comp_time_par] = MSCORparallel(objFun, x0MatSelected,...
    MaxTime, MaxRuns, MaxIter, 1, 2, TolFun1, TolFun2, phi, lambda,...
    1, 2, 1, 0);
Time_taken = toc;

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
M_opt = Patient_Trans_Mat_robust(SeqCell, X, BetaVec_opt, N);

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
Plot_PatientTransitionProb_7_trts(MFinalCell, XRaw, PatientNumbers)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta Coefficients  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BetaIndex : 2 = age_at_diag, 3 = disease_duration, 4 = gender",
%             5 = white, 6 = black.

BetaCellFinal = BetaCell_opt;
for BetaIndex = 2:6
    Plot_Beta_7_trts(BetaIndex,BetaCellFinal,NonNullPositions)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Most frequent treatment transition analyses  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TopHowMany = 16;
[CoeffSummary, SeSummary] = TopHowManyTransFormula_7_trts(TopHowMany, SeqCell, BetaCellFinal,...
    BetaALLBootStdErr, N);

filename = 'Output_CoeffSummary.csv';
writematrix(CoeffSummary, filename);

filename = 'Output_SeSummary.csv';
writematrix(SeSummary, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Age and disease duration effect for Gender X Race %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 = BcD, 2 = GA, 3 = IB, 4 = DF, 5 = Nat, 6 = S1P, 7 = AL
% IB -> S1P (3->6), S1P-> BcD (6->1), IB->DF (3->4), IB->Nat (3->5), 
% DF-> BcD (4->1), DF->S1P (4->6), Nat->BcD (5->1), Nat->S1P (5->6), 
% S1P->DF (6->4).
Numgrids = 100;
Trt_trans = [3 6;
    6 1;
    3 4;
    3 5;
    4 1;
    4 6;
    5 1;
    5 6;
    6 4];

row_labels = {'IB to S1P', 'S1P to BcD', 'IB to DF', 'IB to Nat', 'DF to BcD', ...
              'DF to S1P', 'Nat to BcD', 'Nat to S1P', 'S1P to DF'};
row_labels_initials = {'BcD', 'GA', 'IB', 'DF', 'Nat', 'S1P', 'AL'};
var_names = ["Transition", strcat("grid", string(1:Numgrids))];

%%% Age

min_age = min(XRaw(:,2));
max_age = max(XRaw(:,2));
Ages = linspace(min_age,max_age,Numgrids);
XDummyPatient = nan(Numgrids,p+1);
XDummyPatient(:,1) = 1;
XDummyPatient(:,2) = (Ages - age_mu)./age_sigma;
XDummyPatient(:,3) = (mean(XRaw(:,3))- dd_mu)/dd_sigma;


% (M, W)
XDummyPatient(:,4) = 0; XDummyPatient(:,5) = 1; XDummyPatient(:,6) = 0;
M_dummy_age_MW = DummyPatientsTrans_robust(XDummyPatient, SeqCell, X, BetaCell_opt, N);

T = extract_dynamic_effect(M_dummy_age_MW, N, var_names, row_labels_initials, Trt_trans, row_labels);
writetable(T, 'Effect_age_MW.csv');

% (M, B)
XDummyPatient(:,4) = 0; XDummyPatient(:,5) = 0; XDummyPatient(:,6) = 1;
M_dummy_age_MB = DummyPatientsTrans_robust(XDummyPatient, SeqCell, X, BetaCell_opt, N);

T = extract_dynamic_effect(M_dummy_age_MB, N, var_names, row_labels_initials, Trt_trans, row_labels);
writetable(T, 'Effect_age_MB.csv');

% (F, W)
XDummyPatient(:,4) = 1; XDummyPatient(:,5) = 1; XDummyPatient(:,6) = 0;
M_dummy_age_FW = DummyPatientsTrans_robust(XDummyPatient, SeqCell, X, BetaCell_opt, N);

T = extract_dynamic_effect(M_dummy_age_FW, N, var_names, row_labels_initials, Trt_trans, row_labels);
writetable(T, 'Effect_age_FW.csv');

% (F, B)
XDummyPatient(:,4) = 1; XDummyPatient(:,5) = 0; XDummyPatient(:,6) = 1;
M_dummy_age_FB = DummyPatientsTrans_robust(XDummyPatient, SeqCell, X, BetaCell_opt, N);

T = extract_dynamic_effect(M_dummy_age_FB, N, var_names, row_labels_initials, Trt_trans, row_labels);
writetable(T, 'Effect_age_FB.csv');

%%% Dis_dur

min_dd = min(XRaw(:,3));
max_dd = max(XRaw(:,3));
Numgrids= 100;
DDs = linspace(min_dd,max_dd,Numgrids);
XDummyPatient = nan(Numgrids,p+1);
XDummyPatient(:,1) = 1;
XDummyPatient(:,2) = (mean(XRaw(:,2))- age_mu)/age_sigma;
XDummyPatient(:,3) = (DDs - dd_mu)./dd_sigma;

% (M, W)
XDummyPatient(:,4) = 0; XDummyPatient(:,5) = 1; XDummyPatient(:,6) = 0;
M_dummy_DD_MW = DummyPatientsTrans_robust(XDummyPatient, SeqCell, X, BetaCell_opt, N);

T = extract_dynamic_effect(M_dummy_DD_MW, N, var_names, row_labels_initials, Trt_trans, row_labels);
writetable(T, 'Effect_DD_MW.csv');

% (M, B)
XDummyPatient(:,4) = 0; XDummyPatient(:,5) = 0; XDummyPatient(:,6) = 1;
M_dummy_DD_MB = DummyPatientsTrans_robust(XDummyPatient, SeqCell, X, BetaCell_opt, N);

T = extract_dynamic_effect(M_dummy_DD_MB, N, var_names, row_labels_initials, Trt_trans, row_labels);
writetable(T, 'Effect_DD_MB.csv');

% (F, W)
XDummyPatient(:,4) = 1; XDummyPatient(:,5) = 1; XDummyPatient(:,6) = 0;
M_dummy_DD_FW = DummyPatientsTrans_robust(XDummyPatient, SeqCell, X, BetaCell_opt, N);

T = extract_dynamic_effect(M_dummy_DD_FW, N, var_names, row_labels_initials, Trt_trans, row_labels);
writetable(T, 'Effect_DD_FW.csv');

% (F, B)
XDummyPatient(:,4) = 1; XDummyPatient(:,5) = 0; XDummyPatient(:,6) = 1;
M_dummy_DD_FB = DummyPatientsTrans_robust(XDummyPatient, SeqCell, X, BetaCell_opt, N);

T = extract_dynamic_effect(M_dummy_DD_FB, N, var_names, row_labels_initials, Trt_trans, row_labels);
writetable(T, 'Effect_DD_FB.csv');



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
DummySequences = CreateDummySequence_robust(RandSeed, XDummyPatient, ...
    NumPatientsEachClass, EachSeqLength, SeqCell, X, BetaCell_opt, N);

%%% 
Plot_DummySequences_20Treatments_7_trts(DummySequences,XDummyPatient_Raw)

Time_taken