% Generates an operator Curve

%for use during cmd line run in linux
%addpath '../lib'

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%effective bounds
alpha_effective_bound = alpha / ( 1 - beta);
beta_effective_bound = beta / ( 1 - alpha);

%tuneing parameters
p_prime = 0.9;
delta = 0.03;

%boundaries of accept and reject region
p_one = p_prime + delta;
p_zero = p_prime - delta;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%convinent names for probablity ratios
p_upper = (1 - p_one) / (1 - p_zero);
p_lower = (p_one) / (p_zero);

h = -500:0.01:500;

%operator curve first coordinate (parametric) - p
p = (ones(size(h)) - p_upper.^h) ./ (p_lower.^h - p_upper.^h);

%Graphical test intercepts and slope
h_zero = (log(B)) / (log(p_lower) - log(p_upper));
h_one = (log(A)) / (log(p_lower) - log(p_upper));
s = (log(1 / p_upper)) / (log(p_lower) - log(p_upper));

%Graphical test Paramters and Boundaries.  Purely for grpahing examples, not used in simulation.
m = 0:1:1000;
L_1 = (s .* m) + (h_one .* ones(size(m)));
L_0 = (s .* m) + (h_zero .* ones(size(m)));

%A sample run
set(0,'DefaultAxesFontSize',20)
set(gca, 'FontSize', 20)
f = figure(1);
p_actual = 0.75;
run_sample = count_ones(Sample1dim(p_actual,length(m)));
hold on;
plot(m,L_0,'b--','LineWidth',4);
plot(m,L_1,'b-','LineWidth',4);
plot(m,run_sample,'m-o','LineWidth',4);
xlim([0,50])
ylim([-5,55])
titlestr = sprintf('Sample run of right side test with p_n = %2.2f',p_actual);
xlabel('Sample index (m_n)');
title(titlestr);
legendInfo{1} = 'Lower Boundary L_0';
legendInfo{2} = 'Upper Boundary L_1';
legendInfo{3} = 'Random Walk d_n';
legend(legendInfo);

set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);
