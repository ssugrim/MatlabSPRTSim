% FIXME

%for use during cmd line run in linux
%addpath '../lib'

%%%%%%%%ERROR RATE GLOBALS %%%%%%%
%bounds on alpha effective and beta effective
alpha = 0.05;
beta = 0.01;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%global "time / sample" index, the further out this goes the further out we
%check for crossings.
m = 0:1:100;


%%%%%%%%%% right side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_2 = 0.7;
delta_2 = 0.05;

%boundaries of accept and reject region
p_2_1 = p_prime_2 + delta_2;
p_2_0 = p_prime_2 - delta_2;

%convinent names for probablity ratios
p_upper_2 = (1 - p_2_1) / (1 - p_2_0);
p_lower_2 = (p_2_1) / (p_2_0);

%Graphical test intercepts and slope
h_2_0 = (log(B)) / (log(p_lower_2) - log(p_upper_2));
h_2_1 = (log(A)) / (log(p_lower_2) - log(p_upper_2));
s_2 = (log(1 / p_upper_2)) / (log(p_lower_2) - log(p_upper_2));

%Graphical test Paramters and Boundaries.  Purely for grpahing examples, not used in simulation.
L_2_1 = (s_2 .* m) + (h_one_2 .* ones(size(m)));
L_2_0 = (s_2 .* m) + (h_2_0 .* ones(size(m)));

%%%%%%%%%% left side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_1 = 0.3;
delta_1 = 0.05;

%boundaries of accept and reject region
p_1_1 = p_prime_ls + delta_1;
p_1_0 = p_prime_ls - delta_1;

%convinent names for probablity ratios
p_upper_1 = (1 - p_1_1) / (1 - p_1_0);
p_lower_1 = (p_1_1) / (p_1_0);

%Graphical test intercepts and slope
h_1_0 = (log(B)) / (log(p_lower_1) - log(p_upper_1));
h_1_1 = (log(A)) / (log(p_lower_1) - log(p_upper_1));
s_1 = (log(1 / p_upper_1)) / (log(p_lower_1) - log(p_upper_1));

%Graphical test Paramters and Boundaries.  Purely for grpahing examples, not used in simulation.
L_1_1 = (s_1 .* m) + (h_1_1 .* ones(size(m)));
L_1_0 = (s_1 .* m) + (h_1_0 .* ones(size(m)));

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
        ls_res = SeqThreshold(L_1_1,run_sample,L_1_0);
        rs_res = SeqThreshold(L_2_1,run_sample,L_2_0);
        
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