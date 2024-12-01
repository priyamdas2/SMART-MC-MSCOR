function [x0_opt_sa, fval_sa, comp_time] = saMultSphere(fun_asVec,x0,M,seed)
rng(seed)
tic;
B = size(x0,1);
x0_vec = nan(1,B*M);
for b = 1:B
    x0_vec(((b-1)*M + 1):(b*M)) = x0(b,:)';
end
lb = -1*ones(1,B*M);
ub = 1*ones(1,B*M);
options = optimoptions('simulannealbnd', 'Display', 'off');
[x_opt_vec, fval_sa] = simulannealbnd(fun_asVec,x0_vec,lb,ub, options);

comp_time = toc;

x0_opt_sa = nan(B,M);
for b = 1:B
    x0_opt_sa(b,:) = x_opt_vec(((b-1)*M+1):(b*M))/norm(x_opt_vec(((b-1)*M+1):(b*M)));
end

fprintf('====================== sa starts =============================\n')
fprintf('\n');
fprintf('function value at optimal solution: %d\n', fval_sa);
fprintf('xxxxxxxxxxxxxxxxxxxxxx sa ends xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n')
end

