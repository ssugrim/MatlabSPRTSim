%Computes p_miss, p_fa as a function of p'

%needed when running on linux.
addpath '../lib'

%bounds on alpha effective and beta effective 
%the correct choice of alpha is alpha < beta since our hypothesis
%is not the usual one. alpha will govern or p_miss_effective. 
% When the width of the interval becomes very small, this makes a big
% diffrence.
alpha = 0.01;
beta = 0.05;

%Granularity, smaller takes longer.
gran = 0.01;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%tuneing parameters - Center and width
delta = 0.05;

%have to start a little before the smallest delta, other wise we get odd
%values from the oc curve. 
p_prime = 0.06:.01:0.5;

%storage vars
p_miss = zeros(size(p_prime));
p_fa = zeros(size(p_prime));

for i = 1:1:length(p_prime)
    p_one = p_prime(i) - delta
    
    [p_1,L_1,E_n_1] =  GenOCTheoCurve(alpha, beta, p_prime(i),delta,gran);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %build function to represent L_p
    
    %drop any nan pairs;
    couple = [p_1',L_1'];
    res = couple(0 == sum(isnan(couple), 2), :);
    p_1_c = res(:,1);
    L_1_c = res(:,2);
    
    %Generate a smoothing fit.
    L_1_fobj = fit(p_1_c, L_1_c,'smoothingspline');
    
    %L_2 as a function of p
    L_2_fun = @(p) feval(L_1_fobj,p)+ feval(L_1_fobj,1 - p);
    
    p_fa_eff = integral(L_2_fun, p_one, 1 - p_one,'ArrayValued',true) * (1 / (1 - (2 * p_one)));
    
    p_fa(i) = p_fa_eff;
    %1 - L_1 as a function of p
    L_2_rfun = @(p) 1 - L_2_fun(p);
        
    %effective miss as a numerical integral
    p_miss_eff = integral(L_2_rfun, 0, p_one,'ArrayValued',true) * (2 / p_one)
    p_miss(i) = p_miss_eff;

end

save('./MissVprim.mat');

