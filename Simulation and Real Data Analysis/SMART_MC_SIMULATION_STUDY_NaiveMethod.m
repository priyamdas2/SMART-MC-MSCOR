clear all
DoBootstrapNow = 1; % 1 = Bootstrap is performed now; 0 = Bootstrap s.e. are loaded from saved results

X = readmatrix('Dummy_X.csv');
TrtSeqMat = readmatrix('Dummy_sequences.csv');
NumSubjects = size(X, 1);
p = size(X, 2) - 1;
N = 10;
InterceptIndic = 1;
load('BetaCellTrue.mat', 'BetaCellTrue')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding Bootstrap Standard errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumBootRep = 100; % Dont change; saved result is available only for value 100

if(DoBootstrapNow == 1)
    BetaALLBootStdErr_Naive = BetaALL_BootStdErr_Naive (NumBootRep, X, TrtSeqMat, N);
    filename = ['BetaALL_SIMULATION_BootStdErr_Naive_rep_',num2str(NumBootRep),'.mat'];
    save(filename, 'BetaALLBootStdErr_Naive');
else
    filename = ['BetaALL_SIMULATION_BootStdErr_Naive_rep_',num2str(NumBootRep),'.mat'];
    TempMat = load(filename);
    BetaALLBootStdErr_Naive = TempMat.BetaALLBootStdErr_Naive;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Naive estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RandSeed = 2; % Keeping same as for SMART-MC
rng(RandSeed)
SeqCell = Trt_Seq_Mat2Cell(TrtSeqMat, N);
x0Mat = Generate_Rand_Initial_BetaMat(RandSeed, N, p);
x0Vec = BetaMat2Vec(x0Mat);

objFunNaive = @(BetaVec) -log_likelihood_NaiveMethod(SeqCell, X, BetaVec, N);

len = length(x0Vec); 
lb = -ones(len,1); 
ub = ones(len,1);
x0 = x0Vec; 
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[x_opt, fval_opt] = fmincon(objFunNaive, x0, [], [], [], [], lb, ub, [], options);
BetaVec_opt = x_opt;
BetaCell_opt = Beta_vec2cell(BetaVec_opt, N);

%%% Normalizing each block to have L2 norm 1
BetaCell_opt_scaled = BetaCell_opt;  % Initialize scaled version

[nRows, nCols] = size(BetaCell_opt);  % Get cell array size

for i = 1:nRows
    for j = 1:nCols
        vec = BetaCell_opt{i, j};  % Extract vector
        norm_val = norm(vec, 2);  % Compute L2 norm
        
        if norm_val > 0
            BetaCell_opt_scaled{i, j} = vec / norm_val;  % Normalize
        else
            BetaCell_opt_scaled{i, j} = vec;  % Keep zero vector unchanged
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%

BetaCellEst_Naive = BetaCell_opt_scaled;


percentage = 10;
TopHowManyInitials = round((percentage/100)*N);
TopHowManyTrans = round((percentage/100)*(N^2));
[CoeffsEst_Naive, CoeffsTrue, DiffMatTrans, RMSE, MAD]=Simulation_TopHowManyBetaComparison(TopHowManyInitials,...
    TopHowManyTrans, SeqCell,BetaCellEst_Naive,BetaCellTrue,InterceptIndic,N);
    

filename = 'Output_Simulation_CoeffsEst_Naive.csv';
writematrix(round(CoeffsEst_Naive,2), filename);


NumContVars = 5;
NumIndicVars = 0;
InterceptIndic = 1;

TopTransitions = CoeffsEst_Naive(:,1:2);
Se = zeros(TopHowManyInitials+TopHowManyTrans, NumContVars + 1);

for ii = 1:(TopHowManyInitials+TopHowManyTrans)
    if(ii <= TopHowManyInitials)
        coeffs_here = BetaALLBootStdErr_Naive(1,TopTransitions(ii,1), 1);
        for coeffIdx = 2:(NumContVars + 1)
            coeffs_here = [coeffs_here, BetaALLBootStdErr_Naive(1,TopTransitions(1,1), coeffIdx)];
        end
        Se(ii,:) = coeffs_here;
    else
        coeffs_here = BetaALLBootStdErr_Naive(TopTransitions(ii,1) + 1,TopTransitions(ii,2), coeffIdx);
        for coeffIdx = 2:(NumContVars + 1)
            coeffs_here = [coeffs_here, BetaALLBootStdErr_Naive(TopTransitions(ii,1) + 1,TopTransitions(ii,2), coeffIdx)];
        end
        Se(ii,:) = coeffs_here;
    end
end

filename = 'Output_Simulation_StdErr_Naive.csv';
writematrix(round(Se,3), filename);