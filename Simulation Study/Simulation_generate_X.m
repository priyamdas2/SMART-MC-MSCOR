function X = Simulation_generate_X(RandSeed,NumSubjects,NumContVars,NumIndicVars,BernProb,InterceptIndic,dist)

rng(RandSeed)
X = ones(NumSubjects, 1+NumContVars+NumIndicVars);

if(InterceptIndic == 0)
    X(:,1) = zeros(NumSubjects,1);
end

if(dist == 1)
    X(:,2:(NumContVars+1)) = 20*rand(NumSubjects,NumContVars) - 10;
else
    if(dist == 2)
        X(:,2:(NumContVars+1)) = normrnd(0,1,NumSubjects,NumContVars);
    end
end

if(NumIndicVars > 1)
    X(:,(NumContVars+2):(NumContVars+NumIndicVars+1)) = ...
    rand(NumSubjects, NumIndicVars) < BernProb;
end
end