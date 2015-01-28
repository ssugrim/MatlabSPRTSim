% Plots theoretical OC curve and E(n) curve for varying delta.
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
delta_range = (0.05:0.025:0.15);
split = 0.5;


%% figure var generation
f1 = figure(1);
set(findall(f1,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);
f2 = figure(2);
set(findall(f2,'type','text'),'fontSize',20,'fontWeight','bold','LineWidth',4);

p = 0:0.01:1;
legen = zeros(1,length(delta_range));
plotStyle = {'b+-','ro-','gs-','md-','k^-','c*-'};
%% computed params

for delta_index = 1:length(delta_range)
    
  %  data_index = find(delta == delta_range, 1, 'first');
    delta = delta_range(delta_index);
    sprintf('On delta %f',delta)
    legendInfo{delta_index} = strcat('\delta',sprintf('= %2.2f%',delta));
    
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
    
    figure(f1);
    hold on;
    plot(p,L_2_fun(p),plotStyle{mod(delta_index,length(plotStyle))+1},'LineWidth',2);
    hold off;
    
    
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
    figure(f2);
    hold on;
    plot(p,E_n_2_fun(p),plotStyle{mod(delta_index,length(plotStyle))+1},'LineWidth',2);
    hold off;
end

figure(f1);
legend(legendInfo);
xlabel('Parameter p');
ylabel('L_{2}(p): Two-test probablity of declaring H_0 given p');

figure(f2);
legend(legendInfo);
xlabel('Parameter p');
ylabel('Two-test expected sequence Length');


%matlabpool close;
%fname = strcat('OCTheoCurveDeltaRange-',datestr(now, 'mmmdd.yyyy-HH.MM'),'.mat');
%save(fname);
