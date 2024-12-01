clear all
B = 5; % number of blocks, value can be {1,...,10}U{100}  
M = 5;                      % dimension of the BETA vector
NumExp = 100;
ActivateMSCORParallel = 0; % 1 = parallel, else, NOT parallel
values = [2,3,6,7]; 
for gggg = 1:4 
    which_function = values(gggg);
    
    % Choice of function or user can also define their own function
    % 1 = Negative log of product
    % 2 = Griewank
    % 3 = Negative Sum Squares
    % 4 = Exponential
    % 5 = Easom
    % 6 = Rastrigin
    % 7 = Ackley.
    
    
    %%% Starting point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for rep = 1:NumExp
        fprintf('\n');
        fprintf('Performing experiment number: %d \n',rep);
        fprintf('\n');
        rng(rep)
        starting_point_temp = nan(B,M);
        x0 = nan(B,M);
        for b = 1:B
            starting_point_temp(b,:) = -1 + 2*rand(1,M);
            x0(b,:) = starting_point_temp(b,:)/norm(starting_point_temp(b,:));
            if(b == 1)
                x0_vec =  x0(b,:);
            else
                x0_vec = [x0_vec,x0(b,:)];
            end
        end
        %%%%%% BENCHMARK OBJECTIVE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Negative log of product of absolute values
        if(which_function == 1)
            fun = @(X) Neg_log_prod_of_abs_vals_Sum(X);
            fun_asVec = @(X_vec) Neg_log_prod_of_abs_vals_Sum_asVec(X_vec,M);
            soln_each = (1/sqrt(M))*ones(M,1);
            temp = soln_each';
            soln = repmat(temp,B,1);
        end
        
        % Griewank
        if(which_function == 2)
            fun = @(X) Griewank_Sum(X);
            fun_asVec = @(X_vec) Griewank_Sum_asVec(X_vec,M);
            soln_each = (1/sqrt(M))*ones(M,1);
            temp = soln_each';
            soln = repmat(temp,B,1);
        end
        
        
        % Negative Sum Squares Function
        if(which_function == 3)
            fun = @(X) Neg_sum_square_Sum(X);
            fun_asVec = @(X_vec) Neg_sum_square_Sum_asVec(X_vec,M);
            soln_each = [zeros(M-1,1);1];
            temp = soln_each';
            soln = repmat(temp,B,1);
        end
        
        % Exponential Function
        if(which_function == 4)
            upto_M = M-1;
            fun = @(X) Exponential_Sum(X);
            fun_asVec = @(X_vec) Exponential_Sum_asVec(X_vec,M);
            soln_each = [zeros(M-1,1);1];
            temp = soln_each';
            soln = repmat(temp,B,1);
        end
        
        % Easom
        if(which_function == 5)
            fun = @(X) Easom_Sum(X);
            fun_asVec = @(X_vec) Easom_Sum_asVec(X_vec,M);
            upto_M = 2;
            soln_each = [(1/sqrt(upto_M))*ones(upto_M,1);zeros((M-upto_M),1)];
            temp = soln_each';
            soln = repmat(temp,B,1);
        end
        
        % Rastrigin
        if(which_function == 6)
            fun = @(X) Rastrigin_Sum(X);
            fun_asVec = @(X_vec) Rastrigin_Sum_asVec(X_vec,M);
            soln_each = -(1/sqrt(M))*ones(M,1);
            temp = soln_each';
            soln = repmat(temp,B,1);
        end
        
        % Ackley
        if(which_function == 7)
            fun = @(X) Ackley_Sum(X);
            fun_asVec = @(X_vec)Ackley_Sum_asVec(X_vec,M);
            soln_each = (1/sqrt(M))*ones(M,1);
            temp = soln_each';
            soln = repmat(temp,B,1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        %%% MSCOR
        if(ActivateMSCORParallel == 1)
            [x_opt, fval, comp_time_MSCOR] = MSCORparallel(fun, x0);
        else
            [x_opt, fval, comp_time_MSCOR] = MSCOR(fun, x0);
        end
        
        %%% GA
        seed = rep;
        [x0_opt_ga, fval_ga, comp_time_ga] = gaMultSphere(fun_asVec,x0,M,seed);
        
        %%% SA
        [x0_opt_sa, fval_sa, comp_time_sa] = saMultSphere(fun_asVec,x0,M,seed);
        
        %%% IP
        which_algo = 1; % 1 = interior-point, 2 = sqp, 3 = active-set
        [x0_opt_fmincon_1, fval_fmincon_1, comp_time_1] = fminconMultSphere(fun_asVec,x0,M,which_algo);
        
        %%% SQP
        which_algo = 2; % 1 = interior-point, 2 = sqp, 3 = active-set
        [x0_opt_fmincon_2, fval_fmincon_2, comp_time_2] = fminconMultSphere(fun_asVec,x0,M,which_algo);
        
        %%% Active-set
        which_algo = 3; % 1 = interior-point, 2 = sqp, 3 = active-set
        [x0_opt_fmincon_3, fval_fmincon_3, comp_time_3] = fminconMultSphere(fun_asVec,x0,M,which_algo);
              
        
        f_vals = [fval,fval_ga,fval_sa,fval_fmincon_1,fval_fmincon_2,fval_fmincon_3];
        Times = [comp_time_MSCOR,comp_time_ga,comp_time_sa,comp_time_1, comp_time_2, comp_time_3];
        
        if(rep == 1)
            Summary_fun_vals = f_vals;
            Summary_times = Times;
        else
            Summary_fun_vals = [Summary_fun_vals;f_vals];
            Summary_times = [Summary_times;Times];
        end
    end
    
    
    filename = sprintf('Summary_funVals_%d_B_%d_M_%d_NumExp_%d.csv', which_function,B,M,NumExp);
    writematrix(Summary_fun_vals, filename);
    filename2 = sprintf('Summary_times_%d_B_%d_M_%d_NumExp_%d.csv', which_function,B,M,NumExp);
    writematrix(Summary_times, filename2);
    
    
end

