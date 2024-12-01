function [x0_opt_ga, fval_ga, comp_time] = gaMultSphere(fun_asVec,x0,M,seed)
rng(seed)
tic;
B = size(x0,1);
x0_vec = nan(1,B*M);
for b = 1:B
    x0_vec(((b-1)*M + 1):(b*M)) = x0(b,:)';
end

options = optimoptions('ga', 'Display', 'off');
nonlcon = @(x_vec) sphereConstraints(x_vec,M);
%[x_opt_vec, fval_ga] = ga(fun_asVec, B*M);
[x_opt_vec, fval_ga] = ga(fun_asVec, B*M, [], [], [], [], [], [], nonlcon, options);

comp_time = toc;

x0_opt_ga = nan(B,M);
for b = 1:B
    x0_opt_ga(b,:) = x_opt_vec(((b-1)*M+1):(b*M));
end

fprintf('====================== ga starts =============================\n')
fprintf('\n');
fprintf('function value at optimal solution: %d\n', fval_ga);
fprintf('xxxxxxxxxxxxxxxxxxxxxx ga ends xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n')
end


function [c, ceq] = sphereConstraints(x0_vec,M)
B = length(x0_vec)/M;
c = []; 
if(B == 1)
    ceq = norm(x0_vec(1:M)) - 1;
end
if(B == 2)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1];
end
if(B == 3)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1;...
        norm(x0_vec((2*M+1):(3*M))) - 1];
end
if(B == 4)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1;...
        norm(x0_vec((2*M+1):(3*M))) - 1; norm(x0_vec((3*M+1):(4*M))) - 1];
end
if(B == 5)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1;...
        norm(x0_vec((2*M+1):(3*M))) - 1; norm(x0_vec((3*M+1):(4*M))) - 1;...
        norm(x0_vec((4*M+1):(5*M))) - 1];
end
if(B == 6)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1;...
        norm(x0_vec((2*M+1):(3*M))) - 1; norm(x0_vec((3*M+1):(4*M))) - 1;...
        norm(x0_vec((4*M+1):(5*M))) - 1; norm(x0_vec((5*M+1):(6*M))) - 1];
end
if(B == 7)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1;...
        norm(x0_vec((2*M+1):(3*M))) - 1; norm(x0_vec((3*M+1):(4*M))) - 1;...
        norm(x0_vec((4*M+1):(5*M))) - 1; norm(x0_vec((5*M+1):(6*M))) - 1;...
        norm(x0_vec((6*M+1):(7*M))) - 1];
end
if(B == 8)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1;...
        norm(x0_vec((2*M+1):(3*M))) - 1; norm(x0_vec((3*M+1):(4*M))) - 1;...
        norm(x0_vec((4*M+1):(5*M))) - 1; norm(x0_vec((5*M+1):(6*M))) - 1;...
        norm(x0_vec((6*M+1):(7*M))) - 1; norm(x0_vec((7*M+1):(8*M))) - 1];
end
if(B == 9)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1;...
        norm(x0_vec((2*M+1):(3*M))) - 1; norm(x0_vec((3*M+1):(4*M))) - 1;...
        norm(x0_vec((4*M+1):(5*M))) - 1; norm(x0_vec((5*M+1):(6*M))) - 1;...
        norm(x0_vec((6*M+1):(7*M))) - 1; norm(x0_vec((7*M+1):(8*M))) - 1;...
        norm(x0_vec((8*M+1):(9*M))) - 1];
end
if(B == 10)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1;...
        norm(x0_vec((2*M+1):(3*M))) - 1; norm(x0_vec((3*M+1):(4*M))) - 1;...
        norm(x0_vec((4*M+1):(5*M))) - 1; norm(x0_vec((5*M+1):(6*M))) - 1;...
        norm(x0_vec((6*M+1):(7*M))) - 1; norm(x0_vec((7*M+1):(8*M))) - 1;...
        norm(x0_vec((8*M+1):(9*M))) - 1; norm(x0_vec((9*M+1):(10*M))) - 1];
end

if(B == 100)
    ceq = [norm(x0_vec(1:M)) - 1;norm(x0_vec((M+1):(2*M))) - 1;...
        norm(x0_vec((2*M+1):(3*M))) - 1; norm(x0_vec((3*M+1):(4*M))) - 1;...
        norm(x0_vec((4*M+1):(5*M))) - 1; norm(x0_vec((5*M+1):(6*M))) - 1;...
        norm(x0_vec((6*M+1):(7*M))) - 1; norm(x0_vec((7*M+1):(8*M))) - 1;...
        norm(x0_vec((8*M+1):(9*M))) - 1; norm(x0_vec((9*M+1):(10*M))) - 1;...
        norm(x0_vec((10*M+1):(11*M))) - 1; norm(x0_vec((11*M+1):(12*M))) - 1;...
        norm(x0_vec((12*M+1):(13*M))) - 1; norm(x0_vec((13*M+1):(14*M))) - 1;...
        norm(x0_vec((14*M+1):(15*M))) - 1; norm(x0_vec((15*M+1):(16*M))) - 1;...
        norm(x0_vec((16*M+1):(17*M))) - 1; norm(x0_vec((17*M+1):(18*M))) - 1;...
        norm(x0_vec((18*M+1):(19*M))) - 1; norm(x0_vec((19*M+1):(20*M))) - 1;...
        norm(x0_vec((20*M+1):(21*M))) - 1; norm(x0_vec((21*M+1):(22*M))) - 1;...
        norm(x0_vec((22*M+1):(23*M))) - 1; norm(x0_vec((23*M+1):(24*M))) - 1;...
        norm(x0_vec((24*M+1):(25*M))) - 1; norm(x0_vec((25*M+1):(26*M))) - 1;...
        norm(x0_vec((26*M+1):(27*M))) - 1; norm(x0_vec((27*M+1):(28*M))) - 1;...
        norm(x0_vec((28*M+1):(29*M))) - 1; norm(x0_vec((29*M+1):(30*M))) - 1;...
        norm(x0_vec((30*M+1):(31*M))) - 1; norm(x0_vec((31*M+1):(32*M))) - 1;...
        norm(x0_vec((32*M+1):(33*M))) - 1; norm(x0_vec((33*M+1):(34*M))) - 1;...
        norm(x0_vec((34*M+1):(35*M))) - 1; norm(x0_vec((35*M+1):(36*M))) - 1;...
        norm(x0_vec((36*M+1):(37*M))) - 1; norm(x0_vec((37*M+1):(38*M))) - 1;...
        norm(x0_vec((38*M+1):(39*M))) - 1; norm(x0_vec((39*M+1):(40*M))) - 1;...
        norm(x0_vec((40*M+1):(41*M))) - 1; norm(x0_vec((41*M+1):(42*M))) - 1;...
        norm(x0_vec((42*M+1):(43*M))) - 1; norm(x0_vec((43*M+1):(44*M))) - 1;...
        norm(x0_vec((44*M+1):(45*M))) - 1; norm(x0_vec((45*M+1):(46*M))) - 1;...
        norm(x0_vec((46*M+1):(47*M))) - 1; norm(x0_vec((47*M+1):(48*M))) - 1;...
        norm(x0_vec((48*M+1):(49*M))) - 1; norm(x0_vec((49*M+1):(50*M))) - 1;...
        norm(x0_vec((50*M+1):(51*M))) - 1; norm(x0_vec((51*M+1):(52*M))) - 1;...
        norm(x0_vec((52*M+1):(53*M))) - 1; norm(x0_vec((53*M+1):(54*M))) - 1;...
        norm(x0_vec((54*M+1):(55*M))) - 1; norm(x0_vec((55*M+1):(56*M))) - 1;...
        norm(x0_vec((56*M+1):(57*M))) - 1; norm(x0_vec((57*M+1):(58*M))) - 1;...
        norm(x0_vec((58*M+1):(59*M))) - 1; norm(x0_vec((59*M+1):(60*M))) - 1;...
        norm(x0_vec((60*M+1):(61*M))) - 1; norm(x0_vec((61*M+1):(62*M))) - 1;...
        norm(x0_vec((62*M+1):(63*M))) - 1; norm(x0_vec((63*M+1):(64*M))) - 1;...
        norm(x0_vec((64*M+1):(65*M))) - 1; norm(x0_vec((65*M+1):(66*M))) - 1;...
        norm(x0_vec((66*M+1):(67*M))) - 1; norm(x0_vec((67*M+1):(68*M))) - 1;...
        norm(x0_vec((68*M+1):(69*M))) - 1; norm(x0_vec((69*M+1):(70*M))) - 1;...
        norm(x0_vec((70*M+1):(71*M))) - 1; norm(x0_vec((71*M+1):(72*M))) - 1;...
        norm(x0_vec((72*M+1):(73*M))) - 1; norm(x0_vec((73*M+1):(74*M))) - 1;...
        norm(x0_vec((74*M+1):(75*M))) - 1; norm(x0_vec((75*M+1):(76*M))) - 1;...
        norm(x0_vec((76*M+1):(77*M))) - 1; norm(x0_vec((77*M+1):(78*M))) - 1;...
        norm(x0_vec((78*M+1):(79*M))) - 1; norm(x0_vec((79*M+1):(80*M))) - 1;...
        norm(x0_vec((80*M+1):(81*M))) - 1; norm(x0_vec((81*M+1):(82*M))) - 1;...
        norm(x0_vec((82*M+1):(83*M))) - 1; norm(x0_vec((83*M+1):(84*M))) - 1;...
        norm(x0_vec((84*M+1):(85*M))) - 1; norm(x0_vec((85*M+1):(86*M))) - 1;...
        norm(x0_vec((86*M+1):(87*M))) - 1; norm(x0_vec((87*M+1):(88*M))) - 1;...
        norm(x0_vec((88*M+1):(89*M))) - 1; norm(x0_vec((89*M+1):(90*M))) - 1;...
        norm(x0_vec((90*M+1):(91*M))) - 1; norm(x0_vec((91*M+1):(92*M))) - 1;...
        norm(x0_vec((92*M+1):(93*M))) - 1; norm(x0_vec((93*M+1):(94*M))) - 1;...
        norm(x0_vec((94*M+1):(95*M))) - 1; norm(x0_vec((95*M+1):(96*M))) - 1;...
        norm(x0_vec((96*M+1):(97*M))) - 1; norm(x0_vec((97*M+1):(98*M))) - 1;...
        norm(x0_vec((98*M+1):(99*M))) - 1; norm(x0_vec((99*M+1):(100*M))) - 1];
end



end
