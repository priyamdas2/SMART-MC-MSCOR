B = 5;
n = 5;   % n_b = n
ActivateMSCORParallel = 0; % 1 = parallel, else, NOT parallel

%% Starting point

starting_point_temp = nan(B,n);
x0 = nan(B,n);
for b = 1:B
    starting_point_temp(b,:) = rand(1,n);
    x0(b,:) = starting_point_temp(b,:)/norm(starting_point_temp(b,:));
end
%% Modified Ackley's function

fun = @(X) Ackley_Sum(X);
soln_each = (1/sqrt(n))*ones(n,1);
temp = soln_each';
soln = repmat(temp,B,1); % true global minima

%% MSCOR
if(ActivateMSCORParallel == 1)
    [x_opt, fval, comp_time_MSCOR] = MSCORparallel(fun, x0);
else
    [x_opt, fval, comp_time_MSCOR] = MSCOR(fun, x0);
end

