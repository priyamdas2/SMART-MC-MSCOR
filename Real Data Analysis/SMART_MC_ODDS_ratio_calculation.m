% RUN 'Odds_ration_plot.R' after this to plot them
N = 7;
NonNullPositions = csvread('Output_NonNullPositions.csv');
CHat = csvread('Output_CHat.csv');
BetaVec_opt = csvread('Output_BetaVec_opt.csv');
age_muSigma_dd_muSigma = csvread('Output_age_muSigma_dd_muSigma.csv');
MS_covariates = csvread('Real Data/MS covariates.csv');
TreatmentNamesShort = ["BcD", "GA", "IB","DF","Nat","S1P","AL"];

age_mu = age_muSigma_dd_muSigma(1);
age_sigma = age_muSigma_dd_muSigma(2);
dd_mu = age_muSigma_dd_muSigma(3);
dd_sigma = age_muSigma_dd_muSigma(4);

%% Treatment transition location
% 1 = BcD, 2 = GA, 3 = IB, 4 = DF, 5 = Nat, 6 = S1P, 7 = AL
% IB-> S1P (3->6), S1P-> BcD (6->1), IB->DF (3->4), IB->Nat (3->5), 
% DF-> BcD (4->1), DF->S1P (4->6), Nat->BcD (5->1), Nat->S1P (5->6), 
% S1P->DF (6->4).
Trt_trans = [3 6;
    6 1;
    3 4;
    3 5;
    4 1;
    4 6;
    5 1;
    5 6;
    6 4];
Num_imp_trans = size(Trt_trans,1);

%% Age array
age_arr = [30, 60];
age_arr_scaled = (age_arr - age_mu) ./ age_sigma;

%% Disease duration quartiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dd_raw = MS_covariates(:,3);
quantile_vec = [0.25,0.5,0.75];
num_quantiles = length(quantile_vec);
dd_quartiles_raw = quantile(dd_raw, quantile_vec);
writematrix(dd_quartiles_raw, 'Output_dd_quartiles_raw.csv');
dd_quartiles = (dd_quartiles_raw - dd_mu) ./ dd_sigma;


%% Sample Patients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_patients = nan(12,6,num_quantiles);
OR = nan(7,12,num_quantiles);
OR_with_names = cell(7,12+2,num_quantiles);
for i = 1:num_quantiles
    
    dd_here = dd_quartiles(i);
    
    % Age = 30, Sex = F, Race = White
    sample_patients(1,:,i) = [1, age_arr_scaled(1), dd_here, 1, 1, 0];
    
    % Age = 30, Sex = F, Race = Black
    sample_patients(2,:,i) = [1, age_arr_scaled(1), dd_here, 1, 0, 1];
    
    % Age = 30, Sex = F, Race = Others
    sample_patients(3,:,i) = [1, age_arr_scaled(1), dd_here, 1, 0, 0];
    
    % Age = 30, Sex = M, Race = White
    sample_patients(4,:,i) = [1, age_arr_scaled(1), dd_here, 0, 1, 0];
    
    % Age = 30, Sex = M, Race = Black
    sample_patients(5,:,i) = [1, age_arr_scaled(1), dd_here, 0, 0, 1];
    
    % Age = 30, Sex = M, Race = Others
    sample_patients(6,:,i) = [1, age_arr_scaled(1), dd_here, 0, 0, 0];
    
    % Age = 60, Sex = F, Race = White
    sample_patients(7,:,i) = [1, age_arr_scaled(2), dd_here, 1, 1, 0];
    
    % Age = 60, Sex = F, Race = Black
    sample_patients(8,:,i) = [1, age_arr_scaled(2), dd_here, 1, 0, 1];
    
    % Age = 60, Sex = F, Race = Others
    sample_patients(9,:,i) = [1, age_arr_scaled(2), dd_here, 1, 0, 0];
    
    % Age = 60, Sex = M, Race = White
    sample_patients(10,:,i) = [1, age_arr_scaled(2), dd_here, 0, 1, 0];
    
    % Age = 60, Sex = M, Race = Black
    sample_patients(11,:,i) = [1, age_arr_scaled(2), dd_here, 0, 0, 1];
    
    % Age = 60, Sex = M, Race = Others
    sample_patients(12,:,i) = [1, age_arr_scaled(2), dd_here, 0, 0, 0];
    
    M = Patient_Trans_Mat_FAST(sample_patients(:,:,i), BetaVec_opt, N, NonNullPositions, CHat);
    
    for j = 1:12
        patient_trans_mat = M(2:end,:,j);
        for k = 1:Num_imp_trans
            p_baseline = patient_trans_mat(Trt_trans(k,1), Trt_trans(k,1));
            p_move = patient_trans_mat(Trt_trans(k,1), Trt_trans(k,2));
            OR(k,j,i) = odds_ratio_2to1(p_baseline, p_move);
        end
    end
end

for i = 1:num_quantiles
    AA = [TreatmentNamesShort(Trt_trans(:,1))',...
        TreatmentNamesShort(Trt_trans(:,2))', num2cell(OR(:,:,i))];
    if i==1
        writematrix(AA, 'Output_OR_q25.csv');
    end
    if i==2
        writematrix(AA, 'Output_OR_q50.csv');
    end
    if i==3
        writematrix(AA, 'Output_OR_q75.csv');
    end  
end
% RUN 'Odds_ration_plot.R' after this to plot them