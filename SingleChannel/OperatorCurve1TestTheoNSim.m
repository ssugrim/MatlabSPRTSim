% Generates an operator Curve, both the theoretical and Simulated. These
% can be compared in a plot. Also generates a theoretical and Simulated
% E(n)

%needed when running on linux.
%addpath '../lib'

matlabpool open

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%effective bounds
alpha_effective_bound = alpha / ( 1 - beta);
beta_effective_bound = beta / ( 1 - alpha);

%tuneing parameters - Center and width
p_prime = 0.2;
delta = 0.1;

%boundaries of accept and reject region
p_one = p_prime + delta;
p_zero = p_prime - delta;

%dummy parameter
h = -500:.001:500;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%Theoretical operator curve second coordinate (parametric) - L(p)
L = ((A.^h) - ones(size(h))) ./ (A.^h - B.^h);

%convinent names for probablity ratios
p_upper = (1 - p_one) / (1 - p_zero);
p_lower = (p_one) / (p_zero);

%Theoretical operator curve first coordinate (parametric) - p
p = (ones(size(h)) - p_upper.^h) ./ (p_lower.^h - p_upper.^h);

%Theoretical Expected sequence length (E_p(n))
E_n = ((L .* log(B)) + ((ones(size(h)) - L) .* log(A))) ./ ((p .* log(p_lower)) + ((ones(size(h)) - p) .* log(p_upper)));

%Graphical test intercepts and slope
h_zero = (log(B)) / (log(p_lower) - log(p_upper));
h_one = (log(A)) / (log(p_lower) - log(p_upper));
s = (log(1 / p_upper)) / (log(p_lower) - log(p_upper));

%Graphical test Paramters and Boundaries.  Purely for grpahing examples, not used in simulation.
m = 0:1:1000;
length(m);
L_1 = (s .* m) + (h_one .* ones(size(m)));
L_0 = (s .* m) + (h_zero .* ones(size(m)));

%Simulation paramters
p_test = 0:.01:1;
deadline = 1000000;
num_steps = 1000;

%Simulation Storage vars
calc_E_n = zeros(size(p_test));
calc_p_h1 = zeros(size(p_test));

%for each value of p_test, run simulation to compute the effective p_0 and
%E_p[n]
for j = 1:1:length(p_test);
    
    %dummy data vars;
    steps = zeros(num_steps,1);
    decision = zeros(size(steps));
    
    parfor i = 1:1:length(steps)
        %function that actually generates the data
        cur_thres = GenOCthreshold(p_test(j),s,h_zero,h_one,deadline);
        
        %result analysis / storage
        steps(i)=cur_thres(4);
        if cur_thres(2) > cur_thres(1)
           decision(i)=1;
        end
    end
    
    %status update
    p_test(j)
    
    %per p value store the result
    calc_E_n(j)=mean(steps);
    calc_p_h1(j)=mean(decision);
end

matlabpool close
%generic save
fname = strcat('oc_curve-',datestr(now, 'mmmdd.yyyy-HH.MM'),'.mat');
save(fname);



set(0,'DefaultAxesFontSize',20)
set(gca, 'FontSize', 20)

view_m = 0.5;

%Plot L(P)
figure(1);
plot(p(p < view_m), L(p < view_m),'k--','LineWidth',4);
hold;
plot(p_test(p_test < view_m),1 - calc_p_h1(p_test < view_m),'r-','LineWidth',4);
xlabel('Parameter (p)');
ylabel('L(p)=P(H_0|p)');
title('Operator Characteristic Curve, Truncated to show detail')
legendInfo{1} = sprintf('Theoretical O.C. Curve');
legendInfo{2} = sprintf('Simulated O.C. Curve');
legend(legendInfo);
line([p_prime,p_prime],[0,1])
line([p_one,p_one],[0,1])
line([p_zero,p_zero],[0,1])
xlim([0,view_m]);

%plot E_N(p)
figure(2);
plot(p, E_n,'k--','LineWidth',4);
hold;
plot(p_test,calc_E_n,'r-','LineWidth',2);
xlabel('Parameter (p_n)');
ylabel('E_{p_n}(m)');
title('Expected stopping time Curve')
legendInfo{1} = sprintf('Theoretical Expected Stopping Time');
legendInfo{2} = sprintf('Simulated Expected Stopping Time');
legend(legendInfo);


