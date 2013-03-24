% Generates a sample run for the 2 test comparison. Will draw labeled plot.

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

%global index, this is how far out the test will run
m = 0:1:1000;


%%%%%%%%%% right side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters
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

%tuneing parameters
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


%A sample run
p_actual = 0.55;
dat = Sample1dim(p_actual,length(m));
run_sample = count_ones(dat);
axes();
hold on;

%plotting TODO add axes labels
plot(m,L_1_ls,'--','color','r');
plot(m,L_0_ls,'--','color','g');
plot(m,L_1_rs,'--','color','m');
plot(m,L_0_rs,'--','color','b');
plot(m,run_sample,'.-','color','k');
xlim([0,50])
ylim([0,50])
act_str = sprintf('Sample Run - P_{actual} = %0.5f', p_actual);
Left_str = sprintf('Left Side L_1 P_prime = %0.5f', p_prime_ls);
Right_str = sprintf('Right Side L_1 P_prime = %0.5f', p_prime_rs);
legend(Left_str,'Left Side L_0',Right_str, 'Right Side L_0', act_str)
hold off;