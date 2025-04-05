function [x_opt, fval, comp_time] = MSCOR(objFun, x0, MaxTime, MaxRuns, MaxIter, sInitial, rho, TolFun1, TolFun2, phi,lambda, DisplayUpdate, DisplayEvery, PrintStepSize, PrintSolution)
% MSCOR: Multiple Spherically Constrained Optimization Routine - a BlackBox
%        optimization algorithm!
%
% Inputs:
%   objFun   - Handle to the objective function, takes B x M matrix
%   x0       - Initial guess for the optimization, B x M matrix, each
%              row lies on (M-1) dimensional unit sphere.
%   MaxTime  - maximum allowed execution time in seconds
%   MaxRuns  - number of maximum runs allowed
%   MaxIter  - number of maximum iterations allowed within each run
%   sInitial - initial step-size
%   rho      - step decay rate
%   TolFun1  - tolerance value 1 (controls step decay rate)
%   TolFun2  - tolerance value 1 (controls run restart decision)
%   phi      - minimum allowed step-size
%   lambda   - sparsity threshold
%   DisplayUpdate - 1 = shows iteration update, 0 = doesn't show.
%   DisplayEvery  - Every that many seconds update is displayed only if
%                   DisplayUpdate is set to 1
%   PrintStepSize - 1 = Displays current step-size while updating, 0 = doesn't show. 
%   PrintSolution - 1 = Displays final solution where function is minimized, 0 = doesn't show. 

% Outputs:
%   x_opt   - The optimal solution (minimizer)
%   fval    - Function value at the optimal solution

%%%%%%%% parameters %%%%%%%%%%%%



if (nargin < 4)
    if (nargin < 3)
        MaxTime = 3600;
        if (nargin < 2)
            fprintf('==> Error: Please provide an initial point. If not sure, try out any B x M matrix of random entries. It will be normalized automatically.  \n')
            return;
        end
    end
    MaxRuns = 1000;
    MaxIter = 10000;
    sInitial = 1;
    rho = 2;
    TolFun1 = 10^(-6);
    TolFun2 = 10^(-20);
    phi = 10^(-20);
    lambda = 10^(-6);
    DisplayUpdate = 1;
    DisplayEvery = 2;
    PrintStepSize = 1;
    PrintSolution = 0;
end

fprintf('========================= MSCOR Starts =======================\n')
tic;
B = size(x0,1);
M = size(x0,2);

row_wise_norm = vecnorm(x0, 2, 2);

fprintf('\n')
fprintf('=> Provided intial point will be normalized each row-wise (i.e., unit-sphere-wise), if not normalized already. \n');
if(min(row_wise_norm) == 0)
    fprintf('\n')
    fprintf('==> ERROR: Norm of atleast one unit-spherical block is 0, try different starting point. \n');
    fprintf('\n')
    return;
end

x0_normalized = x0 ./ row_wise_norm;

ValAtInitialPoint  = objFun(x0_normalized);

fprintf('\n');
fprintf('=> Hmm... obj. fun. value at provided initial point is: %d. \n', ValAtInitialPoint);
fprintf('\n');
fprintf('=> Lets BLACKBOX it!!! \n');
fprintf('\n');

adjustment_factor = 10^(20);
Theta_array = zeros(B,M,MaxRuns);
Loop_solution = zeros(MaxRuns, 1);
array_of_values = zeros(MaxIter,1);
last_toc = 0;
break_now = 0;

for iii = 1:MaxRuns
    epsilon = sInitial;
    epsilon_decreasing_factor = rho;
    if(iii == 1)
        Theta = x0_normalized;
    else
        Theta = Theta_array(:,:,iii-1);
    end
    M = size(Theta,2);
    
    
    for i = 1:MaxIter
        if(toc > MaxTime)
            break_now = 1;
            fprintf('=> As requested, MSCOR has been terminated after %.2f seconds :( \n', MaxTime);
            fprintf('\n')
            break;
        end
        current_lh = objFun(Theta);
        
        %%%% Time display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toc_now = toc;
        if(DisplayUpdate == 1)
            if(toc_now - last_toc > DisplayEvery)
                if(PrintStepSize == 1)
                    fprintf('=> Executing Run: %d, iter: %d, current obj. fun. value: %d, current log10(step-size): %.2f. \n', iii, i, current_lh, log10(epsilon));
                else
                     fprintf('=> Executing Run: %d, iter: %d, current obj. fun. value: %d. \n', iii, i, current_lh);
                end
                last_toc = toc_now;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[iii,i/1000,log10(epsilon), current_lh,norm(theta)]
        
        
        total_lh = zeros(2*M,B);
        matrix_update_at_h = zeros(B,M,2*M);
        
        total_lh_alt = zeros(2*M,B);
        matrix_update_at_h_alt = zeros(B,M,2*M);
        
        for b = 1:B
            theta = Theta(b,:);
            for location_number = 1:(2*M)
                Theta_temp = Theta;
                Theta_temp_2 = Theta;
                change_loc  = ceil(location_number/2);
                possibility = theta;
                possibility(change_loc) = 0;
                significant_positions = find(gt(abs(possibility), lambda*ones(1,M)));
                possibility = zeros(1,M);
                possibility_alt = zeros(1,M);
                epsilon_temp = ((-1)^location_number)*epsilon;
                M_here = length(significant_positions)+1;
                if(M_here >= 2)
                    D = (2*sum(theta(significant_positions)))^2 - 4*(M_here-1)*...
                        (2*theta(change_loc)*epsilon_temp+epsilon_temp^2);
                    while(D <0 && abs(epsilon_temp) > phi)
                        epsilon_temp = epsilon_temp/epsilon_decreasing_factor;
                        D = (2*sum(theta(significant_positions)))^2 - 4*(M_here-1)*...
                            (2*theta(change_loc)*epsilon_temp+epsilon_temp^2);
                    end
                    if(D >= 0)
                        possibility(change_loc) = theta(change_loc) + epsilon_temp;
                        possibility_alt(change_loc) = theta(change_loc) + epsilon_temp;
                        a = M_here-1;
                        bb = 2*sum(theta(significant_positions));
                        c = 2*theta(change_loc)*epsilon_temp+epsilon_temp^2;
                        D = bb^2 - 4*a*c;
                        t_here = (-bb + sqrt(D))/(2*a);
                        t_here_alt = (-bb - sqrt(D))/(2*a);
                        possibility(significant_positions) = theta(significant_positions) + t_here;
                        possibility_alt(significant_positions) = theta(significant_positions) + t_here_alt;
                        Theta_temp(b,:) = possibility;
                        Theta_temp_2(b,:) = possibility_alt;
                        total_lh(location_number,b) = objFun(Theta_temp);
                        total_lh_alt(location_number,b) = objFun(Theta_temp_2);
                    else
                        possibility = theta;
                        possibility_alt = theta;
                        total_lh(location_number,b) = current_lh;
                        total_lh_alt(location_number,b) = current_lh;
                    end
                else
                    possibility(change_loc) = round(theta(change_loc));
                    possibility_alt(change_loc) = round(theta(change_loc));
                    Theta_temp(b,:) = possibility;
                    Theta_temp_2(b,:) = possibility_alt;
                    total_lh(location_number,b) = objFun(Theta_temp);
                    total_lh_alt(location_number,b) = objFun(Theta_temp_2);
                end
                
                matrix_update_at_h(b,:,location_number) = possibility;
                matrix_update_at_h_alt(b,:,location_number) = possibility_alt;
            end
        end
        
        
        [minValue_1, linearIndex_1] = min(total_lh(:));
        [row_1, col_1] = ind2sub(size(total_lh), linearIndex_1);  %[location_number, block_no]
        
        [minValue_2, linearIndex_2] = min(total_lh_alt(:));
        [row_2, col_2] = ind2sub(size(total_lh_alt), linearIndex_2); %[location_number, block_no]
        
        if(minValue_1 < minValue_2)
            if(minValue_1 < current_lh)
                updated_block = matrix_update_at_h(col_1,:,row_1); % [block_no, full vec, location num]
                Theta(col_1,:) = updated_block;
            end
            final_value = min(minValue_1,current_lh);
        else
            if(minValue_2 < current_lh)
                updated_block = matrix_update_at_h_alt(col_2,:,row_2);
                Theta(col_2,:) = updated_block;
            end
            final_value = min(minValue_2,current_lh);
        end
        
        
        
        array_of_values(i) =  min(final_value,current_lh);
        
        if(i > 1)
            if(abs(array_of_values(i) - array_of_values(i-1)) < TolFun1)
                if(epsilon > epsilon_decreasing_factor*phi)
                    epsilon = epsilon/epsilon_decreasing_factor;
                else
                    break
                end
            end
        end
        
    end
    
    Theta_array(:,:,iii) = Theta;
    Loop_solution(iii) = objFun(Theta);
    if(iii > 1)
        old_soln = round(Theta_array(:,:,iii-1)*(adjustment_factor))/adjustment_factor;
        new_soln = round(Theta_array(:,:,iii)*(adjustment_factor))/adjustment_factor;
        if(norm(old_soln - new_soln) < TolFun2)
            break
        end
    end
    
    if(break_now == 1)
        break;
    end
end


x_opt = Theta;
fval = objFun(x_opt);
comp_time = toc;

% Final solution
if(PrintSolution == 1)
    fprintf('\n')
    fprintf('=> Final MSCOR solution is: \n');
    disp(x_opt)
end
fprintf('\n')
fprintf('=> Obj. fun. value at MSCOR minima: %d \n',fval);
fprintf('\n')
fprintf('=> Total time taken: %.4f secs.\n',comp_time);

fprintf('xxxxxxxxxxxxxxxxxxxxxx MSCOR ends xxxxxxxxxxxxxxxxxxxxxxxxxx\n')
end