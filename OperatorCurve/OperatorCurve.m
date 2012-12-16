% Generates an operator Curve

addpath '../lib'

%Jpeg params
res = 1600;
font = 25;

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%effective bounds
alpha_effective_bound = alpha / ( 1 - beta);
beta_effective_bound = beta / ( 1 - alpha);

%tuneing parameters
p_prime = 0.9;
delta = 0.01;

%boundaries of accept and reject region
p_one = p_prime + delta;
p_zero = p_prime - delta;

%dummy parameter
h = -500:.01:500;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%operator curve second coordinate (parametric) - L(p)
L = ((A.^h) - ones(size(h))) ./ (A.^h - B.^h);

%convinent names for probablity ratios
p_upper = (1 - p_one) / (1 - p_zero);
p_lower = (p_one) / (p_zero);

%operator curve first coordinate (parametric) - p
p = (ones(size(h)) - p_upper.^h) ./ (p_lower.^h - p_upper.^h);

%Expected sequence length (E_p(n))
E_n = ((L .* log(B)) + ((ones(size(h)) - L) .* log(A))) ./ ((p .* log(p_lower)) + ((ones(size(h)) - p) .* log(p_upper)));

%Graphical test intercepts and slop
h_zero = (log(B)) / (log(p_lower) - log(p_upper));
h_one = (log(A)) / (log(p_lower) - log(p_upper));
s = (log(1 / p_upper)) / (log(p_lower) - log(p_upper));

%Graphical test Paramters and Boundaries. 
m = 0:1:10000;
L_1 = (s .* m) + (h_one .* ones(size(m)));
L_0 = (s .* m) + (h_zero .* ones(size(m)));

%Graphical Simulation parameters
p_test = 0.75:.01:1;
calc_E_n = zeros(size(p_test));
calc_p_h1 = zeros(size(p_test));

%for each value of p_test, run simulation to compute the effective p_0 and
%E[n];
for j = 1:1:length(p_test);
    steps = zeros(500,1);
    decision = zeros(size(steps));
    for i = 1:1:length(steps)
        cur_samples = count_ones(Sample1dim(p_test(j), length(m)));
        cur_thres = SeqThreshold(L_1, cur_samples ,L_0);
        steps(i)=cur_thres(4);
        if cur_thres(2) > cur_thres(1)
           decision(i)=1;
        end
    end
    calc_E_n(j)=mean(steps);
    calc_p_h1(j)=mean(decision);
end

save('./oc_curve.mat');

%plot(p(p > .75), L(p > .75));
%hold;
%plot(p_test,1 - calc_p_h1);