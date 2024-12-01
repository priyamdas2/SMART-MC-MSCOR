objFun = @(X) Ackley_Sum(X);

B  = 4;
M = 5;
rng(1)
starting_point_temp = nan(B,M);
x0 = nan(B,M);
for b = 1:B
    starting_point_temp(b,:) = -1 + 2*rand(1,M);
    x0(b,:) = starting_point_temp(b,:)/norm(starting_point_temp(b,:));
end

[x_opt, fval] = MSCOR(objFun, x0, .5);