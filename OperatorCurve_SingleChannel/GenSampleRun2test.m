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
p_prime_rs = 0.8;
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
p_prime_ls = 0.2;
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
p_actual = 0.75;
dat = Sample1dim(p_actual,length(m));
run_sample = count_ones(dat);
axes();
hold on;

set(0,'DefaultAxesFontSize',20)
set(gca, 'FontSize', 20)
f = figure(1);

%plotting TODO add axes labels
plot(m,L_1_ls,'r-','LineWidth',4);
legendInfo{1} = 'Left Upper Boundary';
plot(m,L_1_rs,'b-','LineWidth',4);
legendInfo{2} = 'Right Upper Boundary';
plot(m,L_0_rs,'b--','LineWidth',4);
legendInfo{3} = 'Right Lower Boundary';
plot(m,L_0_ls,'r--','LineWidth',4);
legendInfo{4} = 'Left Lower Boundary';
plot(m,run_sample,'m-o','LineWidth',4);
legendInfo{5} = 'Random Walk d_n';
legend(legendInfo);
xlim([0,50]);
ylim([-5,55]);
xlabel('Sample index (m_n)');
titlestr = sprintf('Sample run of two sided test with p_n = %2.2f',p_actual);
title(titlestr);
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);
hold off;
