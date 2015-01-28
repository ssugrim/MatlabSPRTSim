% Generates an E(n) curve across the [0,1] interval. Can be compared to the
% E(n) curve of OperatorCurve

%for use during cmd line run in linux
%addpath '../lib'

%%%%%%%%ERROR RATE GLOBALS %%%%%%%
%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%effective bounds
alpha_effective_bound = alpha / ( 1 - beta);
beta_effective_bound = beta / ( 1 - alpha);

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%global "time / sample" index, the further out this goes the further out we
%check for crossings.
m = 0:1:100;


%%%%%%%%%% right side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_rs = 0.7;
delta_rs = 0.05;

%boundaries of accept and reject region
p_one_rs = p_prime_rs + delta_rs;
p_zero_rs = p_prime_rs - delta_rs;

%convinent names for probablity ratios
p_upper_rs= (1 - p_one_rs) / (1 - p_zero_rs);
p_lower_rs = (p_one_rs) / (p_zero_rs);

%Graphical test intercepts and slope
h_zero_rs = (log(B)) / (log(p_lower_rs) - log(p_upper_rs));
h_one_rs = (log(A)) / (log(p_lower_rs) - log(p_upper_rs));
s_rs = (log(1 / p_upper_rs)) / (log(p_lower_rs) - log(p_upper_rs));

%Graphical test Paramters and Boundaries.  Purely for grpahing examples, not used in simulation.
L_1_rs = (s_rs .* m) + (h_one_rs .* ones(size(m)));
L_0_rs = (s_rs .* m) + (h_zero_rs .* ones(size(m)));

%%%%%%%%%% left side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_ls = 0.3;
delta_ls = 0.05;

%boundaries of accept and reject region
p_one_ls = p_prime_ls + delta_ls;
p_zero_ls = p_prime_ls - delta_ls;

%convinent names for probablity ratios
p_upper_ls= (1 - p_one_ls) / (1 - p_zero_ls);
p_lower_ls = (p_one_ls) / (p_zero_ls);

%Graphical test intercepts and slope
h_zero_ls = (log(B)) / (log(p_lower_ls) - log(p_upper_ls));
h_one_ls = (log(A)) / (log(p_lower_ls) - log(p_upper_ls));
s_ls = (log(1 / p_upper_ls)) / (log(p_lower_ls) - log(p_upper_ls));

%Graphical test Paramters and Boundaries.  Purely for grpahing examples, not used in simulation.
L_1_ls = (s_ls .* m) + (h_one_ls .* ones(size(m)));
L_0_ls = (s_ls .* m) + (h_zero_ls .* ones(size(m)));

% breadth of probablities we're going to try. 
p_test = 0:0.01:.99;
num_steps = 1000;

%Simulation Storage vars
calc_E_n = zeros(size(p_test));

%here we reduce to the original miss / false alarm defintion when a actually makes a decision.  
for j = 1:1:length(p_test);
    
    %dummy data vars;
    steps = zeros(num_steps,1);
    decision = zeros(size(steps));
    
    for i = 1:1:length(steps);
        %generate a d vector
        dat = Sample1dim(p_test(j),length(m));
        run_sample = count_ones(dat);
        
        %compute the threshods
        ls_res = SeqThreshold(L_1_ls,run_sample,L_0_ls);
        rs_res = SeqThreshold(L_1_rs,run_sample,L_0_rs);
        
        %smalles index is the first crossing
        steps(i) = min(ls_res(4),rs_res(4));
    end
    %current status reporting
    p_test(j)
    
    %store the mean number of steps for indivudal p
    calc_E_n(j)=mean(steps);
end

%store the results for comparison
save('./fast_combine.mat');