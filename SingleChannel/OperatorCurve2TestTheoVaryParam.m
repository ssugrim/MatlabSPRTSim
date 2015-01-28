% Computes P_fa_eff and P_miss_eff as a function of increasing delta,
% Notes it's effect on max(E(n)).

%needed when running on linux.
%addpath '/home/ssugrim/My Documents/Academics/research/ChannelSurfing/Simulation/SampleGen/lib'
%addpath '../lib'

%matlabpool open;

%% Non SPRT Parameters

%Granularity, smaller takes longer.
gran = 0.01;


%% SPRT paramenters

% one of these should be a range.

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%tuneing parameters - Center and width
p_prime = 0.2;
delta_range = (0.125);
split = 0.5;

%% Storage Vars Across Iterations

%should match size of range from previos sections

p_miss_eff = zeros(1,length(delta_range));
p_fa_eff = zeros(1,length(delta_range));
max_E_n =  zeros(1,length(delta_range));




%% computed params
%parfor delta_index = 1:length(delta_range)
for delta_index = 1:length(delta_range)


  %  data_index = find(delta == delta_range, 1, 'first');
    delta = delta_range(delta_index);
    sprintf('On delta %f',delta)
    
    %Upper and lower sequential test thresholds (aproximated)
    A = (1 - beta) / alpha;
    B = beta / ( 1 - alpha);
    
    %set p_0 to be 30 percent below p_prime
    p_zero = p_prime - (delta * split);
    p_one = p_prime + (delta * (1 - split));
    
    
    %Single test OC Curve
    [p_1,L_1,E_n_1] =  GenOCTheoCurve(alpha, beta, p_one,p_zero,gran);
    
    %reverse p_1
    p_1_r = 1 - p_1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %build function to represent L_p
    
    %drop any nan pairs;
    couple = [p_1',L_1'];
    res = couple(0 == sum(isnan(couple), 2), :);
    p_1_c = res(:,1);
    L_1_c = res(:,2);
    
    % %Generate a smoothing fit.
    L_1_fobj = fit(p_1_c, L_1_c,'smoothingspline');
    
    %L_2 as a function of p - The 2 sided Operator curve.
    L_2_fun = @(p) feval(L_1_fobj,p)+ feval(L_1_fobj,1 - p);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %build a function that represents E_n
    
    %drop any nan pairs;
    couple = [p_1',E_n_1'];
    res = couple(0 == sum(isnan(couple), 2), :);
    p_1_c = res(:,1);
    E_n_1 = res(:,2);
    
    %Generate a smoothing fit.
    E_n_1_fobj = fit(p_1_c, E_n_1,'smoothingspline');
    
    %L_1 as a function of p
    E_n_2_fun = @(p)max(feval(E_n_1_fobj,p),feval(E_n_1_fobj,1 - p));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calculate Miss and False Alarm
    
    %effective false alarm as a numerical integral
    %The H_0 should be p >  1- (p_pmrime_1 - delta) and p < p_pmrime_1 - delta
    
    %false alarm = Type 2 error = declare the null hypothesis is true when it
    %is not = beta
    p_fa_eff(delta_index) = integral(L_2_fun, p_prime, 1 - (p_prime),'ArrayValued',true) * (1 / (1 - (2 * p_prime)));
    
    %1 - L_1 as a function of p
    L_2_rfun = @(p) 1 - L_2_fun(p);
    
    %effective miss as a numerical integral
    p_miss_eff(delta_index) = integral(L_2_rfun, 0, p_prime,'ArrayValued',true) * (2 / p_prime);
    
    p = 0:0.001:1;
    
    max_E_n(delta_index) = max(E_n_2_fun(p));
end

%matlabpool close;
fname = strcat('OCTheoCurveDeltaRange-',datestr(now, 'mmmdd.yyyy-HH.MM'),'.mat');
save(fname);
