% Generates an Theoretical operator Curve for 2-tests. 
%also generates Theoretical E(n)
% Computes P_{fa-effective} and p_{miss-effecitve}

%needed when running on linux.
%addpath '../lib'

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%Granularity, smaller takes longer.
gran = 0.01;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%tuneing parameters - Center and width
p_prime = 0.2;
delta = 0.1;
split = 0.5;

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
p_fa_eff = integral(L_2_fun, p_prime, 1 - (p_prime),'ArrayValued',true) * (1 / (1 - (2 * p_prime)));

%1 - L_1 as a function of p
 L_2_rfun = @(p) 1 - L_2_fun(p);

%effective miss as a numerical integral
p_miss_eff = integral(L_2_rfun, 0, p_prime,'ArrayValued',true) * (2 / p_prime);

%save('./OC2Test.mat');

% plot L(p)
 p = 0:0.001:1;
 f = figure(1);
set(findall(f,'type','text'),'fontSize',14,'fontWeight','bold');
hold on;
plot(p,L_2_fun(p));
xlabel('Paramter p_{actual}');
ylabel('Probablity of test ending in declaration of H_{0}');
titlestr = sprintf('Operator Curve for two simultaneous tests \n p_{fa} = %f \n p_{miss} = %f \n p_1 = %f  \n p_0=%f',p_fa_eff,p_miss_eff,p_one,p_zero);
title(titlestr);
line([p_prime,p_prime],[0,1],'color','r');
hold off;

%Plot E(n)
f = figure(2);
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);
hold on;
plot(p,E_n_2_fun(p));
xlabel('Paramter p_{N}');
ylabel('Stopping Time M_n');
titlestr = sprintf('Expected Number of samples to classifiy');
title(titlestr);
set(findall(f,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);
%line([0,1],[100,100],'color','r','LineWidth',2);
% line([0,1],[50,50],'color','r');
hold off;


%TODO generate the threshold lines for this choice of paramenters.


