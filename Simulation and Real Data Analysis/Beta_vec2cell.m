% Converts Beta row to an (N+1)*N diemnsional cell, with (p+1)
% elements each, filling up ROW-wise across cell blocks
function BetaCell = Beta_vec2cell(BetaVec, N)
p = length(BetaVec)/(N*(N+1)) - 1;
BetaCell = cell(N+1, N);

count = 0;
for i = 1:(N+1)
    for j = 1:N
        BetaCell{i,j} = nan(p+1,1);
        for e = 1:(p+1)
            count = count + 1;
            BetaCell{i,j}(e) = BetaVec(count);
        end
    end
end

